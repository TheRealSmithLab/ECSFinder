

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.*;

import static java.lang.Math.log10;

public class ScanItFast implements Runnable {
    static boolean VERBOSE = false;
    private final String[] key;
    private final ArrayList<String[]> associativeList;
    private static String Path;
    private static String SSZBINARY;
    private double
            shannon = 0;
    private double SSZR_THRESHOLD = -3.0;// alignments scoring below this will be kept (Z-score)
    private int GAP_THRESHOLD = 50;

    private  double sampled_sd;
    private  double sampled_MFE;
    private  double zscore;
    private double sci;

    private double sampled_sd_rc;   // For reverse complement
    private double sampled_MFE_rc;  // For reverse complement
    private double zscore_rc;       // For reverse complement

    private double sci_rc;

    ScanItFast(ArrayList<String[]> associativeList, String[] key, File Path,
               String SSZBINARY, boolean VERBOSE) {
        ScanItFast.Path = Path+"/aln/";
        ScanItFast.SSZBINARY = SSZBINARY;
        ScanItFast.VERBOSE = VERBOSE;
        this.associativeList = associativeList;
        this.key = key;
    }

    public void run() {
        if (VERBOSE)
            System.out.println("- - -> Starting Scan");


        Map<Character, Integer> letterMap = new HashMap<>();
        letterMap.put('A', 0);
        letterMap.put('T', 1);
        letterMap.put('C', 2);
        letterMap.put('G', 3);
        letterMap.put('N', 4);
        letterMap.put('-', 5);
        Map<Character, Integer> letterMapRC = new HashMap<>();
        letterMapRC.put('A', 1);
        letterMapRC.put('T', 0);
        letterMapRC.put('C', 3);
        letterMapRC.put('G', 2);
        letterMapRC.put('N', 4);
        letterMapRC.put('-', 5);
        // remove identical rows or those with too many gaps & N's
        Iterator<String[]> iter = associativeList.iterator();

        ArrayList<int[]> intTab = new ArrayList<>();
        ArrayList<int[]> intTabRC = new ArrayList<>();
        ArrayList<String> UniqueNames = new ArrayList<>();

// Check if the associative list is empty
        if (associativeList.isEmpty()) {
            System.err.println("Error: The associativeList is empty.");
            return;
        }

// Iterate over the associativeList
        while (iter.hasNext()) {
            String[] line = iter.next();
            String sequence = line[1].toUpperCase();
            int[] seqToInt = new int[sequence.length()];
            int[] seqToIntRC = new int[sequence.length()];

            // Convert the sequence characters to integers and their reverse complement
            for (int i = 0; i < sequence.length(); i++) {
                if (letterMap.get(sequence.charAt(i)) == null) {
                    System.err.println("Error: Invalid character in sequence: " + sequence.charAt(i));
                    continue;
                }
                // Convert char to int in the sequence
                seqToInt[i] = letterMap.get(sequence.charAt(i));
                seqToIntRC[sequence.length() - i - 1] = letterMapRC.get(sequence.charAt(i));
            }

            // Split species name and create unique name
            String[] species_part = line[0].split("_");
            String species_part_new = species_part[0] + "_" + species_part[1].split("\\.")[0];

            // If unique name not already present, add the name and corresponding sequences
            if (!UniqueNames.contains(species_part_new)) {
                UniqueNames.add(species_part_new);
                intTabRC.add(seqToIntRC);
                intTab.add(seqToInt);
            }
        }

// Remove rows where 50% or more of the elements are gaps (4) or Ns (5)
        int requiredCount = (int) (0.5 * intTab.get(0).length);

// Create new lists to store filtered rows
        ArrayList<int[]> newArray_gap = new ArrayList<>();
        ArrayList<int[]> newArrayRC_gap = new ArrayList<>();
        ArrayList<String> newUniqueNames = new ArrayList<>();

// Iterate through rows in intTab and intTabRC
        for (int i = 0; i < intTab.size(); i++) {
            int[] originalArray = intTab.get(i);
            int count4 = 0;
            int count5 = 0;

            // Count the number of gaps (4) and Ns (5) in the row
            for (int element : originalArray) {
                if (element == 4) {
                    count4++;
                }
                if (element == 5) {
                    count5++;
                }
            }

            // If fewer than 50% of the elements are gaps/Ns, keep the row and corresponding unique name
            if (count4 < requiredCount && count5 < requiredCount) {
                newArray_gap.add(originalArray);                // Keep row in intTab
                newArrayRC_gap.add(intTabRC.get(i));            // Keep corresponding row in intTabRC
                newUniqueNames.add(UniqueNames.get(i));         // Keep corresponding unique name
            }
        }

// Replace old lists with filtered ones
        intTab.clear();
        intTab.addAll(newArray_gap);

        intTabRC.clear();
        intTabRC.addAll(newArrayRC_gap);

        UniqueNames.clear();
        UniqueNames.addAll(newUniqueNames);

// Remove columns entirely made up of gaps (4) or Ns (5)
        ArrayList<Integer> columnsToRemove = new ArrayList<>();
        int numRows = intTab.size();
        if (numRows<= 3) {
            if (VERBOSE)
                System.out.println("Too many species with gappy sequences");
            return;
        }

        int numCols = intTab.get(0).length;

// Check each column
        for (int col = 0; col < numCols; col++) {
            boolean allGapsOrNs = true;

            for (int row = 0; row < numRows; row++) {
                if (intTab.get(row)[col] != 4 && intTab.get(row)[col] != 5) {
                    allGapsOrNs = false;
                    break;
                }
            }

            if (allGapsOrNs) {
                columnsToRemove.add(col);
            }
        }

// Create new ArrayLists without columns that are entirely made up of gaps/Ns
        ArrayList<int[]> newArray1 = new ArrayList<>();
        ArrayList<int[]> newArray2 = new ArrayList<>();

        for (int[] originalArray : intTab) {
            int[] modifiedArray = new int[numCols - columnsToRemove.size()];
            int colIndex = 0;

            for (int col = 0; col < numCols; col++) {
                if (!columnsToRemove.contains(col)) {
                    modifiedArray[colIndex] = originalArray[col];
                    colIndex++;
                }
            }

            newArray1.add(modifiedArray);
        }

        for (int[] originalArray : intTabRC) {
            int[] modifiedArray = new int[numCols - columnsToRemove.size()];
            int colIndex = 0;

            for (int col = 0; col < numCols; col++) {
                if (!columnsToRemove.contains(col)) {
                    modifiedArray[colIndex] = originalArray[col];
                    colIndex++;
                }
            }

            newArray2.add(modifiedArray);
        }

// Replace old arrays with the filtered ones
        intTab.clear();
        intTab.addAll(newArray1);

        intTabRC.clear();
        intTabRC.addAll(newArray2);

// Convert UniqueNames to an array if needed
        String[] nameTab = new String[UniqueNames.size()];
        nameTab = UniqueNames.toArray(nameTab);

// Check if enough sequences are left
        int goodSeqs = intTab.size();

        if (intTab.size() <= 3) {
            if (VERBOSE)
                System.out.println("-> Not Enough seqs ");
            return;
        }

        if (intTab.size() <= 3) {
            if(VERBOSE)
                System.out.println("-> Not Enough seqs in this window!");
            return;
        }

        if (!(nameTab[0].contains("homo"))) {
            if(VERBOSE)
                System.out.println("-> No human in this alignment block");
            return;
        }



        /*********************************************************************
         calculate stats						*
         *********************************************************************/
        if (VERBOSE)
            System.out.println("- - -> calculating statistics");
        double uniqueSeqs = goodSeqs;
        int outCols = intTab.get(0).length; // Change last variable if CLUSTAL properties change
        double[] stats = new double[6];
        double[] totalChars;
        // boolean[] isNotUnique = new boolean[outCols]; // Update to match the number of columns
        double[] chars;
        double[] column = new double[outCols];

        // Calculate MPI
        for (int k = 0; k < outCols; k++) {
            double identicalNuc = 0.0;
            double totalNuc = 0.0;
            double[][] stats1 = new double[6][6];


            for (int i = 0; i < goodSeqs; i++) {
                for (int j = i + 1; j < goodSeqs; j++) {
                    int charI = intTab.get(i)[k];
                    int charJ = intTab.get(j)[k];

                    if (isValidNucleotide(charI) && isValidNucleotide(charJ)) {
                        stats1[charI][charJ] += 1.0;
                        totalNuc += 1.0;

                        if (charI == charJ) {
                            identicalNuc += 1.0;
                        }
                    }
                }
            }
            if (totalNuc > 0) {
                column[k] = identicalNuc / totalNuc;
            } else {
                column[k] = 0.0; // Avoid division by zero
            }

        }

        double sum = 0.0;
        for (double value : column) {
            sum += value;
        }

        double newMPI = sum / column.length;


        totalChars = new double[]{0.0, 0.0, 0.0, 0.0, 0.0};

        for (int i = 0; i < intTab.get(0).length; i++) {
            chars = new double[]{0.0, 0.0, 0.0, 0.0, 0.0};
            for (int j = 0; j < goodSeqs; j++) {
                if (intTab.get(j)[i] == 5) {
                    chars[4] += 1.0;
                    totalChars[4] += 1.0;
                } else {
                    chars[intTab.get(j)[i]] += 1.0;
                    totalChars[intTab.get(j)[i]] += 1.0;
                }
            }
            for (int z = 0; z != 5; z++) {
                double probz = chars[z] / uniqueSeqs;
                shannon = (chars[z] == 0) ? shannon + 0 : shannon + probz * (Math.log(probz) / Math.log(2));
            }

        }
        Map<Integer, Character> reverseMap = new HashMap<>();
        reverseMap.put(0, 'A');
        reverseMap.put(1, 'T');
        reverseMap.put(2, 'C');
        reverseMap.put(3, 'G');
        reverseMap.put(4, 'N');
        reverseMap.put(5, '-');

        // prepare clustalw file
        if (VERBOSE)
            System.out.println("- -> preparing Clustal format");
        String[] outAln = new String[goodSeqs];
        String[] outAlnRC = new String[goodSeqs];
        int iterate = 0;
        for (int seq = 0; seq < intTab.size(); seq++) { //removed x < goodseqs

            outAln[iterate] = nameTab[seq].substring(0, Math.min(nameTab[seq].length(), 20));
            outAlnRC[iterate] = nameTab[seq].substring(0, Math.min(nameTab[seq].length(), 20));
            for (int i = 0; i != 25 - Math.min(nameTab[seq].length(), 20); i++) {
                outAln[iterate] = outAln[iterate] + " ";
                outAlnRC[iterate] = outAlnRC[iterate] + " ";
            }
            for (int i = 0; i != intTab.get(0).length; i++) {

                outAln[iterate] = outAln[iterate] + reverseMap.get(intTab.get(seq)[i]);

                outAlnRC[iterate] = outAlnRC[iterate] + reverseMap.get(intTabRC.get(seq)[i]);

            }
            outAln[iterate] = outAln[iterate] + "\n";
            outAlnRC[iterate] = outAlnRC[iterate] + "\n";
            iterate++;


        }

        stats[0] = newMPI * 100;                                                                        // Mean Pairwise ID
        double var = 0.0;
        for (double v : column) {
            var = var + Math.pow(v - newMPI, 2);
        }
        double standard = Math.sqrt(var / (column.length - 1));
        stats[1] = standard;                        // Variance
        stats[2] = -1 * shannon / ((double) outCols);       // Normalized Shanon entropy
        stats[3] = 100 * (totalChars[2] + totalChars[3]) / (totalChars[0] + totalChars[1] + totalChars[2] + totalChars[3]);       // GC content
        stats[4] = 100 * totalChars[4] / (outCols * goodSeqs);// GAP content


        // save BED coords
        if (VERBOSE)
            System.out.println("- -> Calculating BED coords ");
        String BedFile = key[0] + "\t";
        BedFile = BedFile + key[1] + "\t" + key[2] + "\t";


        BedFile = BedFile + goodSeqs + "_" + ((double) (int) (10 * stats[0]) / 10) + "_"      // MPI
                + ((double) (int) (10 * stats[1]) / 10) + "_"            // STDEV
                + ((double) (int) (100 * stats[2]) / 100) + "_"                // SHANNON
                + ((double) (int) (10 * stats[3]) / 10) + "_"                  //      GC
                + ((double) (int) (10 * stats[4]) / 10);                     // GAPS


        if (VERBOSE)
            System.out.println("Pre SISSIz bed file: \n" + " " + BedFile);
        int random = (int) ((double) 10000 * Math.random());
        File Aln = new File(Path  + BedFile.replaceAll("\t", "_") + ".aln." + random),    //
                AlnRC = new File(Path  + BedFile.replaceAll("\t", "_") + "rc.aln." + random);  //
        // v v v v v v v v    INCLUSION STATS     v v v v v v v v v v v v v
        //MPI greater or equal than 50 and Gap content smaller or equal than 50
        if (stats[4] <= GAP_THRESHOLD && stats[0] >= 50) {
            // Write Sequences to ALN Format
            try (BufferedWriter WriteClustal = new BufferedWriter(new FileWriter(Aln));
                 BufferedWriter WriteClustalRC = new BufferedWriter(new FileWriter(AlnRC))) {
                WriteClustal.write("CLUSTAL format \n\n");
                WriteClustalRC.write("CLUSTAL format \n\n");

                for (int y = 0; y != goodSeqs; y++) {

                    WriteClustal.write(outAln[y]);
                    WriteClustalRC.write(outAlnRC[y]);
                }
                WriteClustal.close();
                WriteClustalRC.close();
            } catch (IOException Err) {
                if (VERBOSE)
                    System.err.println("Arrgh... Couldn't write clustal file!");
                Err.printStackTrace();
                Aln.delete();
                AlnRC.delete();
                return;
            }



        } else {
            if (VERBOSE) {
                System.out.println("---> rejected alignment");
                System.out.println("     outcols = " + outCols + "\tuniqueseqs = " + uniqueSeqs +
                        "\tGAPS = " + stats[4] + "\n    PID = " + stats[0]);
                if (stats[0] < 5)
                    System.out.println("-----> SUPER LOW PID");
            }
            Aln.delete();
            AlnRC.delete();
            return;
        }
        String FinalBedFile,
                FinalBedFileRC,
                Antisense = (key[3].equals("+")) ? "-" : "+";


        //***************** 	SISSIz scan & parse		******************
        String[] SissizOutTab = new String[12];


        try {

            SissizOutTab = ScanSSZ(Path, BedFile, random);

            if (SissizOutTab == null) { // timeout
                Aln.delete();
            }
        } catch (IOException Err) {
            Err.printStackTrace();
            System.err.println("ScanSSZ failed with ");

            Aln.delete();
        }
        // delete empty alignments
        if (SissizOutTab == null || SissizOutTab[12] == null) {
            Aln.delete();
        } else {
            FinalBedFile = BedFile + "_" + (int) (Double.parseDouble(SissizOutTab[12]) * -100) + "_" + key[3];
            // delete low scoring alignments
            if (Double.parseDouble(SissizOutTab[12]) > SSZR_THRESHOLD) {


                Aln.delete();


            } else {
                sci= Double.parseDouble(SissizOutTab[7]);
                sampled_MFE = Double.parseDouble(SissizOutTab[10]);
                sampled_sd = Double.parseDouble(SissizOutTab[11]);
                zscore = Double.parseDouble(SissizOutTab[12]);
                String fileNameBed = FinalBedFile.replace("\t", "_");
                File theDir = new File(Path + "/" + fileNameBed);
                if (!theDir.exists()) {
                    theDir.mkdirs();
                }


                File NewFile = new File(theDir, fileNameBed+".aln");
                int file_count = 0;
                while (NewFile.exists()) {
                    file_count++;
                    NewFile = new File(theDir,fileNameBed + ".aln_" + file_count);
                }
                boolean result = Aln.renameTo(NewFile);
                // Run RNAalifold
                String clustalFilePath = NewFile.getAbsolutePath();

                runRNAalifold(clustalFilePath,fileNameBed, String.valueOf(theDir));

                boolean rScapeResult = runRScape(fileNameBed + ".stk", String.valueOf(theDir));

// If R-scape returns false, skip alignment
                if (!rScapeResult) {
                    System.err.println("Skipping alignment because R-scape failed or encountered a fatal error.");

                    return;
                }
                FilterOutput filterOutput = new FilterOutput();

                double eval = filterOutput.processFilesWithSuffix(String.valueOf(theDir), "helixcov", "E-value: ");
                double cov= filterOutput.processFilesWithSuffix(String.valueOf(theDir), "power", "# BPAIRS observed to covary ");
                double[] energies = filterOutput.processTxtFiles(String.valueOf(theDir));
                int len_prediction =Integer.valueOf(key[2])-Integer.valueOf(key[1]);
                double[] array_variates= new double[]{
                        energies[0] / len_prediction,
                        energies[1] / len_prediction,
                        log10(eval),
                        cov,
                        ((double) (int) (10 * stats[0]) / 10),
                        sampled_MFE / len_prediction,
                        sampled_sd,
                        zscore,
                        sci,
                };

                try {
                    // Copy .aln files before deleting the directory
                    copyAlnFiles(theDir, Path);
                    deleteDirectory(theDir);
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }

                String path_csv = ECSFinder.OUT_PATH+"/csv/"+fileNameBed+".csv";

                try {
                    writeFeaturesToCSV(array_variates, path_csv,String.valueOf(theDir));
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }

            }
        }
        // * * * * * *  now for the RC  * * * * * *
        try {
            SissizOutTab = ScanSSZ(Path, BedFile + "rc", random);
            if (SissizOutTab == null) {
                AlnRC.delete();
            }
        } catch (IOException Err) {
            Err.printStackTrace();
            System.err.println("ScanSSZ failed in RC with ");

            AlnRC.delete();
        }
        if (SissizOutTab == null || SissizOutTab[12] == null) {
            AlnRC.delete();
        } else {
            FinalBedFileRC = BedFile + "_" + (int) (Double.parseDouble(SissizOutTab[12]) * -100) + "_" + Antisense;
            // delete low scoring alignments
            if (Double.parseDouble(SissizOutTab[12]) > SSZR_THRESHOLD) {
                AlnRC.delete();

                //    System.out.println(FinalBedFileRC.replaceAll("_", "\t"));
            } else {
                sci_rc= Double.parseDouble(SissizOutTab[7]);
                sampled_MFE_rc = Double.parseDouble(SissizOutTab[10]);
                sampled_sd_rc = Double.parseDouble(SissizOutTab[11]);
                zscore_rc = Double.parseDouble(SissizOutTab[12]);
                String fileNameBedRc = FinalBedFileRC.replace("\t", "_");
                // Create the necessary directory for the alignment
                File theDir = new File(Path + "/" + fileNameBedRc);
                if (!theDir.exists()) {
                    theDir.mkdirs();
                }
                File NewFile = new File(theDir,  fileNameBedRc + ".aln");
                int file_count = 0;
                while (NewFile.exists()) {
                    file_count++;
                    NewFile = new File(theDir,  fileNameBedRc +  ".aln_" + file_count);
                }
                boolean result = AlnRC.renameTo(NewFile);
                // Run RNAalifold
                String clustalFilePath = NewFile.getAbsolutePath();



                runRNAalifold(clustalFilePath,fileNameBedRc, String.valueOf(theDir));

                // Run R-scape

                boolean rScapeResult = runRScape(fileNameBedRc + ".stk", String.valueOf(theDir));

// If R-scape returns false, skip alignment
                if (!rScapeResult) {
                    System.err.println("Skipping alignment because R-scape failed or encountered a fatal error.");

                    return;
                }
                FilterOutput filterOutput = new FilterOutput();
                double eval = filterOutput.processFilesWithSuffix(String.valueOf(theDir), "helixcov", "E-value: ");
                double cov= filterOutput.processFilesWithSuffix(String.valueOf(theDir), "power", "# BPAIRS observed to covary ");

                double[] energies = filterOutput.processTxtFiles(String.valueOf(theDir));
                int len_prediction =Integer.valueOf(key[2])-Integer.valueOf(key[1]);
                double[] array_variates= new double[]{

                        energies[0] / len_prediction,
                        energies[1] / len_prediction,
                        log10(eval),
                        cov,
                        ((double) (int) (10 * stats[0]) / 10),
                        sampled_MFE_rc / len_prediction,
                        sampled_sd_rc,
                        zscore_rc,
                        sci_rc

                };

                try {
                    // Copy .aln files before deleting the directory
                    copyAlnFiles(theDir, Path);

                    deleteDirectory(theDir);
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
                String path_csv = ECSFinder.OUT_PATH+"/csv/"+fileNameBedRc+".csv";
                try {
                    writeFeaturesToCSV(array_variates, path_csv, String.valueOf(theDir));
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }

            }
            return;
        }

    }

    /*********************************************************************
     SISSIz scan & parse						*
     //*********************************************************************/
    // sissiz-di       cluster.109999_step.aln  8       150     0.8759  0.8542  0.0094  -13.88  -8.20   3.48    -1.63
    protected static String[] ScanSSZ(String path, String bedFile, int id) throws IOException {
        String name = path + "/" + bedFile.replaceAll("\t", "_") + ".aln." + id;
        List<String> command = Arrays.asList(SSZBINARY, "-j", "-t", "--sci", name);

        long timeoutMs = 300_000;  // 5 minutes
        List<String> output;
        try {
            output = ECSFinder.runExternalCommand(command, new File(path), timeoutMs, VERBOSE);
        } catch (InterruptedException e) {
            e.printStackTrace();
            return null;
        } catch (IOException e) {
            e.printStackTrace();
            return null;
        }

        // parse the output
        String[] sissizOutTab = new String[12];
        for (String line : output) {
            if (line.startsWith("TREE")) {
                String[] parts = line.split(";");
                if (parts.length > 1 && parts[1].startsWith("sissiz")) {
                    if (VERBOSE) {
                        System.out.println(parts[1]);
                    }
                    sissizOutTab = parts[1].split("\\s");
                }
            }
        }
        return sissizOutTab;
    }



    //*********************************************************************
    //						Sample process						*
    //*********************************************************************
    private static boolean isAlive(Process p) {
        try {
            p.exitValue();
            return false;
        } catch (IllegalThreadStateException e) {
            return true;
        }
    }

    public void setSszR(double newValue) {
        SSZR_THRESHOLD = newValue;
    }
    private static String getDirectoryPath(String filePath) {
        Path path = Paths.get(filePath);
        Path parentPath = path.getParent();
        return parentPath != null ? parentPath.toAbsolutePath().toString() : "";
    }
    public void setGap(int newGap) {
        GAP_THRESHOLD = newGap;
    }

    // Helper function to check if a character is a valid nucleotide
    private boolean isValidNucleotide(int c) {
        // valid nucleotides are 0, 1, 2, 3, and 5 (excluding 4)
        return (c >= 0 && c <= 5) && (c != 4);
    }
    // Method to extract the file name without the extension
    private static String getFileNameWithoutExtension(String filePath) {
        Path path = Paths.get(filePath);
        String fileName = path.getFileName().toString();
        return fileName.replaceFirst("[.][^.]+$", "");
    }
    // Function to run RNAalifold
    private void runRNAalifold(String clustalFile, String noExt, String directoryPath) {
        // Build the RNAalifold command
        List<String> cmd = Arrays.asList(
                ECSFinder.RNAALIFOLD,
                "--noLP",
                "-r",
                "--noPS",
                "--aln-stk=" + noExt,
                "--id-prefix=" + noExt,
                clustalFile
        );


        long timeoutMs = 300_000;  // 5 * 60 * 1000

        try {

            List<String> outputLines = ECSFinder.runExternalCommand(
                    cmd,
                    new File(directoryPath),
                    timeoutMs,
                    VERBOSE
            );

            // If you want to capture stdout in a .txt file:
            File outputTextFile = new File(directoryPath, noExt + ".txt");
            try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputTextFile))) {
                for (String line : outputLines) {
                    writer.write(line);
                    writer.newLine();
                }
            }

        } catch (IOException | InterruptedException e) {
            e.printStackTrace();
        }
    }



    /**
     * Run R-scape on the given .stk file inside directoryPath.
     *
     * @param stkFile        The filename of the .stk file (without the directory prefix).
     * @param directoryPath  The directory where the .stk file resides.
     * @return true if R-scape completed without fatal errors, false otherwise.
     */
    private boolean runRScape(String stkFile, String directoryPath) {
        // Build the R-scape command
        List<String> cmd = Arrays.asList(
                "R-scape",
                "--lancaster",
                "--nofigures",
                "-s",
                directoryPath + "/" + stkFile
        );


        long timeoutMs = 300_000;  // 5 minutes

        try {
            // Run the command with the helper method
            List<String> outputLines = ECSFinder.runExternalCommand(
                    cmd,
                    new File(directoryPath),
                    timeoutMs,
                    VERBOSE
            );

            // Check the lines for any fatal error markers
            boolean fatalErrorFound = false;
            for (String line : outputLines) {
                // If R-scape encounters a fatal error, mark it for skipping
                if (line.contains("Fatal exception") || line.contains("IncompleteGamma")) {
                    System.err.println("R-scape encountered a numerical issue and will skip this alignment.");
                    fatalErrorFound = true;
                }
            }

            // If a fatal error was found, or an exception was thrown
            // due to a non-zero exit code,
            // we handle that here:
            return !fatalErrorFound;

        } catch (IOException | InterruptedException e) {
            // If there's an error or a non-zero exit code, we log and return false
            e.printStackTrace();
            return false;
        }
    }




    // Copy .aln files to the target directory
    public static void copyAlnFiles(File sourceDir, String targetDir) throws IOException {
        // Convert the target directory from String to Path
        Path targetPath = Paths.get(targetDir);

        // Ensure the target directory exists, or create it
        if (!Files.exists(targetPath)) {
            Files.createDirectories(targetPath);
        }

        // Get list of .aln files in the source directory
        File[] alnFiles = sourceDir.listFiles((dir, name) -> name.endsWith(".aln"));

        // Copy each .aln file to the target directory
        if (alnFiles != null) {
            for (File alnFile : alnFiles) {
                Path targetFile = targetPath.resolve(alnFile.getName());
                Files.copy(alnFile.toPath(), targetFile, StandardCopyOption.REPLACE_EXISTING);
            }
        }
    }

    public static void deleteDirectory(File directory) throws IOException {
        if (!directory.exists()) {
            return;
        }

        File[] files = directory.listFiles();
        if (files != null) {
            for (File file : files) {
                if (file.isDirectory()) {
                    deleteDirectory(file);
                } else {
                    if (!file.delete()) {
                        throw new IOException("Failed to delete file: " + file);
                    }
                }
            }
        }

        if (!directory.delete()) {
            throw new IOException("Failed to delete directory: " + directory);
        }
    }

    private static void writeFeaturesToCSV(double[] data, String csvPath, String nameFile) throws IOException {
        boolean fileExists = new File(csvPath).exists();
        BufferedWriter writer = new BufferedWriter(new FileWriter(csvPath, true)); // Enable append mode
            String fileName = nameFile.substring(nameFile.lastIndexOf('/') + 1);
            writer.write((fileName) + ",");
            for (int i = 0; i < data.length; i++) {
                writer.write(String.valueOf(data[i]));
                if (i < data.length - 1) {
                    writer.write(",");
                }
            }
            writer.write("\n");
            writer.close();
        }


}

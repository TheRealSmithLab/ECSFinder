
import java.io.*;
import java.util.*;
import java.util.concurrent.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class ECSFinder {

    static int GAPS = 50, NTHREDS = 4;
    static boolean VERBOSE = false;
    static String FILENAME = "", OUT_PATH = "", dirProgram = "";
    static String SSZBINARY = "", ALIFOLDBINARY = "",
            RNAALIFOLD="",R= "", RSCAPE="";
    static double SSZR = -3.0;
    private static File path;

    public static void main(String[] args) throws IOException, InterruptedException {

        // usage info
        if (args.length == 0) {
            System.out.println("\n\t\t\t  version 1.0 \n" +
                    " ________    ______   ______    ________  _                 __                \n" +
                    "|_   __  | .' ___  |.' ____ \\  |_   __  |(_)               |  ]               \n" +
                    "  | |_ \\_|/ .'   \\_|| (___ \\_|   | |_ \\_|__   _ .--.   .--.| | .---.  _ .--.  \n" +
                    "  |  _| _ | |        _.____`.    |  _|  [  | [ `.-. |/ /'`\\' |/ /__\\\\[ `/'`\\] \n" +
                    " _| |__/ |\\ `.___.'\\| \\____) |  _| |_    | |  | | | || \\__/  || \\__., | |     \n" +
                    "|________| `.____ .' \\______.' |_____|  [___][___||__]'.__.;__]'.__.'[___]    \n" +
                    "\t SCAN MULTIPLE ALIGNMENTS FOR CONSERVED RNA STRUCTURES\n\n" +
                    "Reads a set of maf files, calculates stats, scans with SISSIz and R-scape , outputs bed coordinates of high-confidence predictions\n\n" +
                    "Usage: java ECSFinder [options] -o output/directory -i input.maf (last parameter must be -i)\n\n" +
                    "Output: Two types of results are produced:\n" +
                    " (1) the multiple sequence alignments associated with significant predictions \n" +
                    "     are saved to files in the folder specified with the \"-o\" option.\n" +
                    "     File names correspond to their genomic coordinates in a .bed-compatible format. Ex:\n" +
                    "     output/directory/chrX_12345_12500_80:75:23:14:8:z_300_+.aln\n" +
                    " (2) The genomic coordinates (.bed format) of ECSs are also written to the SDOUT\n" +
                    " (see additional options below).\n" +

                    "Options:\n" +
                    " -c int number of CPUs for calculations (default 4)\n" +
                    " -g int max gap percentage of sequences for 2D prediction (default 50)\n" +
                    " -sszr double report SISSIz+RIBOSUM hits below this Z-score (default -3.0)\n" +
                    " -v verbose (messy but detailed) output\n");
            System.exit(0);
        }

        // get binary paths
        setBinaryPaths();

        // parse arguments
        parseArguments(args);

        // preprocess MAF file using MergeNFilter
        preprocessMafFiles();

        // run RNALalifold and process results
        runRNALalifoldAndProcessResults();
        // Now check for the existence of the output file
        String nameOutputECS= "predicted_ECS.csv";

        callRScript(OUT_PATH + "/structure_input_sense.csv", OUT_PATH + "/"+nameOutputECS);
        callRScript(OUT_PATH + "/structure_input_antisense.csv", OUT_PATH + "/"+nameOutputECS);

        File outputFile = new File(OUT_PATH + "/"+nameOutputECS);
        if (!outputFile.exists()) {
            System.out.println("No ECS were found.");
            System.exit(0); // Exit the program if the file doesn't exist
        } else {
            try {
                // Open the CSV file
                BufferedReader reader = new BufferedReader(new FileReader(outputFile));
                String line;

                // Read the header line and skip
                reader.readLine();
                int first=0;
                // Process each line
                while ((line = reader.readLine()) != null) {
                    // Split the line by commas
                    String[] parts = line.split(",");

                    // Check if Predicted_Class is "TP"
                    if (parts.length >= 3 && "TP".equals(parts[2])) {
                        // Replace underscores with tabs in name_file (the first column)
                        String nameFileWithTabs = parts[0].replace('_', '\t');
                        if(first==0){
                            System.out.println("chrm\tstart\tend\tNum_species\tMPI\tsd\tmean_shannon\tgc\tgap\tzscore\tstrand\tprob");
                        }
                        // Print the name_file and probability to the terminal
                        System.out.println(nameFileWithTabs + "\t" + parts[1]);
                        first+=1;
                    }
                    else if (parts.length >= 3 && "FP".equals(parts[2])) {
                        // Use the original file name
                        String nameFileWithUnderscores = parts[0];

                        // Construct the full path to the .aln file
                        File fileToDelete = new File(path, nameFileWithUnderscores + ".aln");

                        // Check if the file exists and delete it
                        if (fileToDelete.exists()) {
                            fileToDelete.delete();
                        }
                    }


                }

                // Close the reader
                reader.close();
            } catch (IOException e) {
                e.printStackTrace();
            }

            // Also delete the structure_input_sense.csv file
            File senseCsvFile = new File(OUT_PATH, "structure_input_sense.csv");
            File antisenseCsvFile = new File(OUT_PATH, "structure_input_antisense.csv");
            // Check if the CSV file exists and delete it
            if (senseCsvFile.exists()) {
                senseCsvFile.delete();
            }
            if (antisenseCsvFile.exists()){
                antisenseCsvFile.delete();
            }

        }

    }



    private static void setBinaryPaths() throws IOException {
        ALIFOLDBINARY = getBinaryPath("RNALalifold");
        SSZBINARY = getBinaryPath("SISSIz");
        RNAALIFOLD = getBinaryPath("RNAalifold");
        R = getBinaryPath("R");
        RSCAPE = getBinaryPath("R-scape");
    }

    private static String getBinaryPath(String binaryName) throws IOException {
        ProcessBuilder pbCmd = new ProcessBuilder("which", binaryName);
        Process pbCmdProcess = pbCmd.start();
        BufferedReader reader = new BufferedReader(new InputStreamReader(pbCmdProcess.getInputStream()));
        String binaryPath = reader.readLine();
        if (binaryPath == null) {
            System.out.println("Please install " + binaryName + " and link it to your $PATH");
            System.exit(0);
        }
        reader.close();
        return binaryPath;
    }

    private static void parseArguments(String[] args) {
        for (int i = 0; i < args.length; i++) {
            switch (args[i]) {
                case "-c":
                    NTHREDS = Integer.parseInt(args[i + 1]);
                    i++;
                    break;
                case "-g":
                    GAPS = Integer.parseInt(args[i + 1]);
                    i++;
                    break;
                case "-o":
                    OUT_PATH = System.getProperty("user.dir") + "/" + args[i + 1];
                    createDirectory(OUT_PATH);
                    dirProgram = System.getProperty("user.dir");
                    i++;
                    break;
                case "-v":
                    VERBOSE = true;
                    break;
                case "-sszr":
                    SSZR = Double.parseDouble(args[i + 1]);
                    i++;
                    break;
                case "-i":
                    FILENAME = args[i + 1];
                    break;
            }
        }
    }

    private static void createDirectory(String path) {
        File dir = new File(path);
        if (!dir.isDirectory()) {
            dir.mkdirs();
        }
    }

    private static void preprocessMafFiles() {
        File dir = new File(FILENAME);
        if (!dir.isDirectory()) {
            System.out.println("The provided path is not the full path directory name.");
            System.exit(0);
        }

        File[] mafFiles = dir.listFiles((d, name) -> name.endsWith(".maf"));
        if (mafFiles == null || mafFiles.length == 0) {
            System.out.println("No .maf files found in the directory.");
            System.exit(0);
        }

        String[] filePaths = Arrays.stream(mafFiles).map(File::getAbsolutePath).toArray(String[]::new);

        try {
            MergeNFilter mergeNFilter = new MergeNFilter();
            mergeNFilter.process(filePaths, OUT_PATH);
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(0);
        }
    }

    private static void runRNALalifoldAndProcessResults() throws IOException, InterruptedException {
        // Extract file path
        String inputFile = OUT_PATH + "/output.maf";
        String[] nameAlifold = inputFile.split("\\.");
        String Path2 = OUT_PATH + "/stockholm" + nameAlifold[nameAlifold.length - 1];
        File stockholmFolder = new File(Path2);
        if (!stockholmFolder.exists()) {
            stockholmFolder.mkdir();
        }
        executeCommand(inputFile, nameAlifold);

        BufferedReader readFile = new BufferedReader(new FileReader(inputFile));
        int blockAln = 0;
        StringBuilder temp = new StringBuilder();

        ExecutorService multiThreads = Executors.newFixedThreadPool(NTHREDS);
        List<Future<?>> futures;

        String line;
        while ((line = readFile.readLine()) != null) {
            if (line.length() >= 1 && line.charAt(0) != '#') {
                if (line.charAt(0) == 'a') {
                    blockAln++;
                } else if (line.charAt(0) == 's') {
                    temp.append(line).append("@");
                }

            } else if ((temp.toString().split("@").length <= 2) && line.equals("")) {
                temp = new StringBuilder();
            } else if ((temp.toString().split("@").length >= 3) && line.equals("")) { // at least 3 sequences
                String[] tempTab = temp.toString().split("@");
                temp = new StringBuilder();

                ArrayList<String[]> associativeList = new ArrayList<>();
                String[] mafTabTemp = tempTab[0].split("\\s+");
                futures = new ArrayList<>();
                path = new File(OUT_PATH + "/aln/");
                if (!path.isDirectory()) {
                    path.mkdirs();  // Create the directory if it doesn't exist
                }

                int numbZeros;
                String gcReference = "";
                String gcSScons = "";
                String lastDigit = "";
                StringBuilder finalName = new StringBuilder();
                lastDigit += String.valueOf(blockAln);
                numbZeros = 4 - lastDigit.length();
                int x = 0;
                while (x < numbZeros) {
                    finalName.append(0);
                    x++;
                }
                finalName.append(lastDigit);

                File file = new File(OUT_PATH + "/stockholm" +
                        nameAlifold[nameAlifold.length - 1] + "/alifold_"
                        + finalName + ".stk");

                if (file.exists()) {
                    try {
                        if (VERBOSE)
                            System.out.println("Stockholm file: " +
                                    file.getAbsolutePath().substring(file.getAbsolutePath().lastIndexOf("/")
                                            + 1));
                        BufferedReader reader = new BufferedReader(new FileReader(file));
                        String currentLine = reader.readLine();

                        String[] arrayName = new String[5];

                        while (currentLine != null) {

                            if (currentLine.startsWith("#=GF ID ")) {
                                arrayName = currentLine.split("[_.]");

                                associativeList = new ArrayList<>();
                            } else if (currentLine.startsWith("#=GC RF")) {
                                String[] lineReference = currentLine.split(" ");
                                gcReference = lineReference[lineReference.length - 1];
                            } else if (currentLine.startsWith("#=GC SS_cons")) {
                                String[] lineReference = currentLine.split(" ");
                                gcSScons = lineReference[lineReference.length - 1];
                            } else if (currentLine.startsWith("#") || currentLine.equals("")) {
                                currentLine = reader.readLine();
                                continue;

                            } else if (!(currentLine.startsWith("//"))) {

                                String[] species = currentLine.split(" ", 2);
                                species[1] = species[1].trim();

                                associativeList.add(species);
                            }
                            if ((!associativeList.isEmpty()) && currentLine.startsWith("//")) {
                                int[] cordMotif = getRealCoordinates(Integer.parseInt(arrayName[3])
                                        , mafTabTemp, associativeList.get(0)[1]);
                                String loci = Arrays.toString(cordMotif);
                                String chrom = mafTabTemp[1].substring((mafTabTemp[1].lastIndexOf(".")) + 1);
                                String lociChrm = chrom + ", " + loci.substring(1, loci.length() - 1) + ", " +
                                        mafTabTemp[4] + ", " + arrayName[3] + ", " + arrayName[4] + ", " + gcReference + ", "
                                        + gcSScons;
                                String[] arrayLociChrm = lociChrm.split(", ");

                                if (Integer.parseInt(arrayLociChrm[2]) - Integer.parseInt(arrayLociChrm[1]) < 50) {
                                    currentLine = reader.readLine();
                                    continue;
                                }
                                ScanItFast aln = new ScanItFast(associativeList, arrayLociChrm, path,
                                        SSZBINARY, VERBOSE);
                                aln.setSszR(SSZR);
                                aln.setGap(GAPS);

                                Future<?> f = multiThreads.submit(aln);
                                futures.add(f);
                            }
                            currentLine = reader.readLine();
                        }

                        reader.close();

                    } catch (FileNotFoundException e) {
                        e.printStackTrace();
                    }

                } else {
                    System.out.println("File: " + file.getAbsolutePath().substring(file.getAbsolutePath().lastIndexOf("/")
                            + 1) + " does not exist");
                }
                for (Future<?> g : futures) {
                    try {
                        g.get();

                    } catch (Exception Err) {
                        System.err.println("MultiThreads took too long!  OOps!");
                        Err.printStackTrace();
                    }
                }

            }

        }
        // Delete file with stockholm info
        String[] entries = stockholmFolder.list();
        assert entries != null;
        if (entries.length > 0) {
            for (String s : entries) {
                File currentFile = new File(stockholmFolder.getPath(), s);
                currentFile.delete();
            }
        }
        stockholmFolder.delete();

        readFile.close();
        multiThreads.shutdown();
        multiThreads.awaitTermination(60 * 10L, TimeUnit.SECONDS);


    }





    private static int[] getRealCoordinates(int start, String[] mafCord, String motifHuman) {
        int[] cordFinal;
        int[] cordFinalPlus1 = new int[2];
        String withoutGap = mafCord[6].substring(0, start);
        int nuc = withoutGap.replaceAll("-", "").length();

        if (mafCord[4].equals("-")) {
            int lociEnd = (Integer.parseInt(mafCord[5]) + 1 - (Integer.parseInt(mafCord[2]) + nuc)) + 1;
            int lociStart = lociEnd - motifHuman.replaceAll("-", "").length();
            cordFinal = new int[]{lociStart, lociEnd};
        } else {
            int lociStart = Integer.parseInt(mafCord[2]) + nuc;
            int lociEnd = lociStart + motifHuman.replaceAll("-", "").length();
            cordFinal = new int[]{lociStart, lociEnd};
        }
        cordFinalPlus1[0] = cordFinal[0];
        cordFinalPlus1[1] = cordFinal[1] - 1;

        return cordFinalPlus1;
    }

    private static void executeCommand(final String file, String[] nameAlifold) {
        String Path = OUT_PATH + "/stockholm" + nameAlifold[nameAlifold.length - 1];

        try {

            ProcessBuilder pb = new ProcessBuilder(ALIFOLDBINARY, "--id-prefix=alifold", "--noLP",
                    "--maxBPspan=300", "--ribosum_scoring", "--aln-stk", file);
            String cmd = ALIFOLDBINARY + " --id-prefix=alifold" + " --noLP" +
                    " --maxBPspan=300" + " --ribosum_scoring" + " --aln-stk " + file;
            if (VERBOSE)
                System.out.println("Executing command " + cmd);

            pb.directory(new File(Path));
            Process process = pb.start();
            InputStream error = process.getInputStream();
            for (int i = 0; i < error.available(); i++) {
                System.out.println("" + error.read());
            }

            BufferedReader reader =
                    new BufferedReader(new InputStreamReader(process.getInputStream()));

            while ((reader.readLine()) != null) {

            }
            process.waitFor();
        } catch (IOException | InterruptedException e) {
            e.printStackTrace();
        }

    }
    public static List<Double> callRScript(String inputCsv, String outputCsv) throws IOException, InterruptedException {

        // Get the path to the R script
        String rScriptPath = getRScriptPath();
        rScriptPath=rScriptPath.replace("target/", "");
        List<String> command = Arrays.asList("Rscript", rScriptPath+"/prediction_RF.R", inputCsv, outputCsv, rScriptPath);

        ProcessBuilder pb = new ProcessBuilder(command);
        Process process = pb.start();

        // Capture standard output and error
        StringBuilder output = new StringBuilder();
        BufferedReader stdOutput = new BufferedReader(new InputStreamReader(process.getInputStream()));
        BufferedReader stdError = new BufferedReader(new InputStreamReader(process.getErrorStream()));

        String s;
        while ((s = stdOutput.readLine()) != null) {
            output.append(s).append("\n");
        }
        while ((s = stdError.readLine()) != null) {
            output.append(s).append("\n");
        }

        process.waitFor();


        // Read the output predictions from the CSV
        BufferedReader reader = new BufferedReader(new FileReader(outputCsv));
        List<Double> predictions = new ArrayList<>();
        String line;
        // Skip the header line
        reader.readLine();
        while ((line = reader.readLine()) != null) {
            String[] parts = line.split(",");
            if (parts.length > 1) {
                try {
                    predictions.add(Double.parseDouble(parts[1].trim()));
                } catch (NumberFormatException e) {
                    e.printStackTrace();
                }
            }
        }
        reader.close();

        return predictions;
    }
    private static String getRScriptPath() throws IOException {
        // Find the directory where the JAR or class is located
        File jarFile = new File(ECSFinder.class.getProtectionDomain().getCodeSource().getLocation().getPath());
        String jarDir = jarFile.getParent(); // Directory of the running JAR or classpath

        // If running from JAR, get the relative path to the script from the JAR directory
        if (jarFile.isFile()) { // Running from a JAR
            return new File(jarDir, "r_scripts/").getAbsolutePath();
        } else { // Running from class, so it's in the project directory
            return new File("r_scripts/").getAbsolutePath();
        }
    }

}


import java.util.*; import java.io.*;
import java.lang.*;

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

    ScanItFast(ArrayList<String[]> associativeList, String[] key, String Path,
               String SSZBINARY, boolean VERBOSE) {
        ScanItFast.Path = Path;
        ScanItFast.SSZBINARY = SSZBINARY;
        ScanItFast.VERBOSE = VERBOSE;
        this.associativeList = associativeList;
        this.key = key;
    }

    public void run() {
        if (VERBOSE)
            System.out.println("- - -> Starting Scan");


        Map <Character, Integer> letterMap = new HashMap<>();
        letterMap.put('A', 0);
        letterMap.put('T', 1);
        letterMap.put('C', 2);
        letterMap.put('G', 3);
        letterMap.put('N', 4);
        letterMap.put('-', 5);
        Map <Character, Integer> letterMapRC = new HashMap<>();
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
        while (iter.hasNext()) {
            String[] line = iter.next();
            String sequence = line[1].toUpperCase();
            int[] seqToInt = new int[sequence.length()];
            int[] seqToIntRC = new int[sequence.length()];


            for (int i = 0; i < sequence.length(); i++) {
                if (letterMap.get(sequence.charAt(i)) ==null){
                    System.out.println("No sequence is found");
                }
                // convert char to int in the sequence
                seqToInt[i] = letterMap.get(sequence.charAt(i));
                seqToIntRC[sequence.length() - i - 1] = letterMapRC.get(sequence.charAt(i));

            }

            String [] species_part = line[0].split("_");
            String species_part_new =species_part[0]+"_"+species_part[1].split("\\.")[0];
            if (! UniqueNames.contains(species_part_new)) {
                UniqueNames.add(species_part_new);
                intTabRC.add(seqToIntRC);
                intTab.add(seqToInt);
            }
        }
// Remove columns entirely made up of 1s or when 50% or more of a row is made up of gaps and Ns
        int requiredCount = (int) (0.5 * intTab.get(0).length);

        // Create a new ArrayList without arrays containing 50% or more 4s
        ArrayList<int[]> newArray_gap = new ArrayList<>();
        ArrayList<int[]> newArrayRC_gap = new ArrayList<>();
        for (int i = 0; i < intTab.size(); i++) {
            int[] originalArray = intTab.get(i);
            int count4 = 0;
            int count5 = 0;

            for (int element : originalArray) {
                if (element == 4) {
                    count4++;
                } if(element ==5){
                    count5++;
                }
            }

            if (count4 < requiredCount && count5 < requiredCount) {
                newArray_gap.add(originalArray);
                newArrayRC_gap.add(intTabRC.get(i));
            }
        }

        intTab.clear();
        intTab.addAll(newArray_gap);
        intTabRC.clear();
        intTabRC.addAll(newArrayRC_gap);
        // Find columns with only the value 1 in either ArrayList
        ArrayList<Integer> columnsToRemove = new ArrayList<>();
        int numRows1 = intTab.size();
        int numRows2 = intTabRC.size();
        int numCols = intTab.get(0).length;

        for (int col = 0; col < numCols; col++) {
            boolean allOnes = true;

            for (int row = 0; row < numRows1; row++) {
                if (intTab.get(row)[col] != 4 ||intTab.get(row)[col] != 5) {
                    allOnes = false;
                    break;
                }
            }

            if (allOnes) {
                columnsToRemove.add(col);
            } else {
                // Check arrayList2 for columns with only the value 1
                for (int row = 0; row < numRows2; row++) {
                    if (intTabRC.get(row)[col] != 4 ||intTabRC.get(row)[col] != 5) {
                        allOnes = false;
                        break;
                    }
                }

                if (allOnes) {
                    columnsToRemove.add(col);
                }
            }
        }

        // Create new ArrayLists without columns containing only 1
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
        intTab.clear();
        intTab.addAll(newArray1);

        intTabRC.clear();
        intTabRC.addAll(newArray2);

        String[] nameTab = new String[UniqueNames.size()];
        nameTab = UniqueNames.toArray(nameTab);
// first check for > 2 seqs
        int goodSeqs = intTab.size();

        if (intTab.size() <= 3) {
            if (VERBOSE)
                System.out.println("-> Not Enough seqs ");
            return;
        }
//check whether theres homo sapiens
     //   if (!(String.valueOf(UniqueNames.get(0)).startsWith("homo_sapiens"))){
     //       return;
     //   }

        // System.out.println("this is the end of block");
        if (intTab.size() <= 3) {
            if (VERBOSE)
                System.out.println("-> Not Enough seqs in this window!");
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


        totalChars = new double[]{0.0,0.0,0.0,0.0,0.0};

        for (int i = 0; i < intTab.get(0).length; i++){
            chars = new double[]{0.0,0.0,0.0,0.0,0.0};
            for(int j = 0; j < goodSeqs; j++){
                if(intTab.get(j)[i] == 5){
                    chars[4]+= 1.0;
                    totalChars[4]+= 1.0;
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
        reverseMap.put(0,'A');
        reverseMap.put(1,'T');
        reverseMap.put(2,'C');
        reverseMap.put(3,'G');
        reverseMap.put(4,'N');
        reverseMap.put(5,'-');

        // prepare clustalw file
        if (VERBOSE)
            System.out.println("- -> preparing Clustal format");
        String[] outAln = new String[goodSeqs];
        String[] outAlnRC = new String[goodSeqs];
        int iterate = 0;
        for (int seq = 0; seq < intTab.size() ; seq++) { //removed x < goodseqs

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
        double standard = Math.sqrt(var/(column.length -1));
        stats[1] = standard;                        // Variance
        stats[2] = -1 * shannon / ((double) outCols);       // Normalized Shanon entropy
        stats[3] = 100 * (totalChars[2] + totalChars[3]) / (totalChars[0] + totalChars[1] + totalChars[2] + totalChars[3]);       // GC content
        stats[4] = 100 * totalChars[4] / (outCols * goodSeqs);// GAP content




        // save BED coords
        if (VERBOSE)
            System.out.println("- -> Calculating BED coords ");
        String BedFile = key[0] + "\t";
        BedFile = BedFile + key[1] + "\t" + key[2] + "\t";

        BedFile = BedFile + goodSeqs + ":" + ((double) (int) (10 * stats[0]) / 10) + ":"      // MPI
                + ((double) (int) (10 * stats[1]) / 10) + ":"            // STDEV
                + ((double) (int) (100 * stats[2]) / 100) + ":"                // SHANNON
                + ((double) (int) (10 * stats[3]) / 10) + ":"                  //      GC
                + ((double) (int) (10 * stats[4]) / 10);                     // GAPS




        if (VERBOSE)
            System.out.println("Pre SISSIz bed file: \n" + " " + BedFile);
        int random = (int) ((double) 10000 * Math.random());
        File Aln = new File(Path + "/" + BedFile.replaceAll("\t", "_") + ".aln." + random),    //
                AlnRC = new File(Path + "/" + BedFile.replaceAll("\t", "_") + "rc.aln." + random);  //
        // v v v v v v v v    INCLUSION STATS     v v v v v v v v v v v v v
        //MPI greater or equal than 50 and Gap content smaller than 75
        if ( stats[4] <= GAP_THRESHOLD  && stats[0] >= 50){
            // Write Sequences to ALN Format
            try {
                BufferedWriter WriteClustal = new BufferedWriter(new FileWriter( Aln )),
                        WriteClustalRC = new BufferedWriter(new FileWriter( AlnRC ));
                WriteClustal.write("CLUSTAL format \n\n") ;
                WriteClustalRC.write("CLUSTAL format \n\n") ;

                for (int y = 0; y != goodSeqs; y++ ) {

                    WriteClustal.write( outAln[ y ] ) ;
                    WriteClustalRC.write( outAlnRC[ y ] ) ;
                }
                WriteClustal.close() ;
                WriteClustalRC.close() ;
            } catch (IOException Err) {
                if (VERBOSE)
                    System.err.println("Arrgh... Couldn't write clustal file!");
                Err.printStackTrace();
                Aln.delete() ;
                AlnRC.delete() ;
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
                Antisense = (key[3].equals("+"))? "-" : "+";


//***************** 	SISSIz scan & parse		******************
        String[] SissizOutTab = new String[12];


        try {
            SissizOutTab = ScanSSZ( Path, BedFile, random);

            if (SissizOutTab == null) { // timeout
                Aln.delete();
            }
        } catch (IOException Err) {
            Err.printStackTrace();
            System.err.println("ScanSSZ failed with ");

            Aln.delete();
        }
        // delete empty alignments
        if (SissizOutTab == null || SissizOutTab[10] == null) {
            Aln.delete();
        } else {
            FinalBedFile = BedFile + "_" + (int) (Double.parseDouble(SissizOutTab[10]) * -100) + "_" + key[3];
            // delete low scoring alignments
            if (Double.parseDouble(SissizOutTab[10]) > SSZR_THRESHOLD) {



                Aln.delete();


            } else {
                //write bed and rename alignment
                System.out.println(FinalBedFile.replaceAll("_", "\t")+"\t"+SissizOutTab[4]+"\t"+
                        SissizOutTab[5]+"\t"+SissizOutTab[6]+"\t"+SissizOutTab[7]+"\t"+SissizOutTab[8]+"\t"+
                        SissizOutTab[9]+"\t"+SissizOutTab[10]+"\t" +key[7]);

                File NewFile = new File(Path + "/" + FinalBedFile.replaceAll("\t", "_") + ".aln");
                int file_count = 0;
                while (NewFile.exists()) {
                    file_count++;
                    NewFile = new File(Path + "/" + FinalBedFile.replaceAll("\t", "_") + ".aln_" + file_count);
                }
                boolean result = Aln.renameTo(NewFile);
            }
        }
        // * * * * * *  now for the RC  * * * * * *
        try {
            SissizOutTab = ScanSSZ( Path, BedFile + "rc", random);
            if (SissizOutTab == null) {
                AlnRC.delete();
            }
        } catch (IOException Err) {
            Err.printStackTrace();
            System.err.println("ScanSSZ failed in RC with ");
            for (int y = 0; y != goodSeqs; y++) {
                //  if (!isNotUnique[y]) {
                //   System.err.println(outAln[y]);
                // }
            }
            AlnRC.delete();
        }
        if (SissizOutTab == null || SissizOutTab[10] == null) {
            AlnRC.delete();
        } else {
            FinalBedFileRC = BedFile + "_" + (int) (Double.parseDouble(SissizOutTab[10]) * -100) + "_" + Antisense;
            // delete low scoring alignments
            if (Double.parseDouble(SissizOutTab[10]) > SSZR_THRESHOLD) {
                AlnRC.delete();

                //    System.out.println(FinalBedFileRC.replaceAll("_", "\t"));
            } else {

                //write bedRC and rename alignment
                System.out.println(FinalBedFileRC.replaceAll("_", "\t")+"\t"+SissizOutTab[4]+"\t"+
                        SissizOutTab[5]+"\t"+SissizOutTab[6]+"\t"+SissizOutTab[7]+"\t"+SissizOutTab[8]+"\t"+
                        SissizOutTab[9]+"\t"+SissizOutTab[10]+"\t" +key[7]);
                File NewFile = new File( Path + "/" + FinalBedFileRC.replaceAll("\t", "_") + ".aln");
                int file_count = 0;
                while (NewFile.exists()) {
                    file_count++;
                    NewFile = new File( Path + "/" + FinalBedFileRC.replaceAll("\t", "_") + ".aln_" + file_count);
                }
                boolean result = AlnRC.renameTo(NewFile);
            }
            return;
        }

    }


    /*********************************************************************
     SISSIz scan & parse						*
     //*********************************************************************/
    // sissiz-di       cluster.109999_step.aln  8       150     0.8759  0.8542  0.0094  -13.88  -8.20   3.48    -1.63
    protected static String[] ScanSSZ (String Path, String BedFile, int id ) throws
            IOException {
        //stats[0] Mean Pairwise ID
        //stats[1] Variance
        //stats[2] Normalized Shannon entropy
        //stats[3] GC content
        //stats[4] GAP content
        String[] SissizOutTab = new String[12];
        String Output, Error = "";
        String Command = SSZBINARY;

        Command = Command + " -j -t " + Path + "/" + BedFile.replaceAll("\t", "_") + ".aln." + id; // RIBOSUM scoring

        try {
            long now = System.currentTimeMillis();
            long timeoutInMillis = 1000L * 300;                          // max 5 minutes
            long finish = now + timeoutInMillis;
            // launch initial SISSIz call
            String name = Path+ "/" + BedFile.replaceAll("\t", "_")+ ".aln."+ id;
            ProcessBuilder pb = new ProcessBuilder(SSZBINARY, "-j", "-t", name);
            Process Sissiz = pb.start();
            BufferedReader SissizErr = new BufferedReader(new InputStreamReader(Sissiz.getErrorStream()));
            if (VERBOSE)
                System.out.println(": Running " + Command);
            while (isAlive(Sissiz)) {
                Thread.sleep(100);
                if (System.currentTimeMillis() > finish) {
                    if (VERBOSE)
                        System.out.println("SISSIz failed to run within time :-(");
                    SissizErr.close();
                    Sissiz.destroy();
                    return null;


                }
            }
            SissizErr.close();
            // get Output if process didn't complete in recursion
            if (SissizOutTab[0] == null) {
                BufferedReader SissizOut = new BufferedReader(new InputStreamReader(Sissiz.getInputStream()));
                while ((Output = SissizOut.readLine()) != null) {
                    // Check if the line starts with "TREE"
                    if (Output.startsWith("TREE")) {
                        String[] Output_new = Output.split(";");
                        if (Output_new.length > 1 && Output_new[1].startsWith("sissiz")) {
                            if (VERBOSE) {
                                System.out.println(Output_new[1]);
                            }
                            SissizOutTab = Output_new[1].split("\\s");
                            // You can modify elements in SissizOutTab as needed
                            // Example: SissizOutTab[1] = "r";
                        }
                    }
                }
                SissizOut.close();
            }


        } catch (Exception err) {
            System.out.println(" Not enough nucleotides in the column " + Command + "\n  counter--> " );
            System.err.println("Not enough nucleotides in the column " + Command + "\n  counter--> " );
            err.printStackTrace();
            System.err.println("===============================");
        }
        return SissizOutTab;
    }

    //*********************************************************************
    //						Sample process						*
    //*********************************************************************
    private static boolean isAlive( Process p ) {
        try {
            p.exitValue();
            return false;
        } catch (IllegalThreadStateException e) {
            return true;
        }
    }
    public void setSszR(double newValue){
        SSZR_THRESHOLD = newValue;
    }
    public void setGap(int newGap){
        GAP_THRESHOLD = newGap;
    }
    // Helper function to check if a character is a valid nucleotide
    private boolean isValidNucleotide(int c) {
        // valid nucleotides are 0, 1, 2, 3, and 5 (excluding 4)
        return (c >= 0 && c <= 5) && (c != 4);
    }


}

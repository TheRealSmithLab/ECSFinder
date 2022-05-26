import java.util.*; import java.io.*;
import java.lang.*;

public class ScanItFast implements Runnable {
    static boolean VERBOSE = false,
            PRINTALL = false,
            RSCAPE = false;
    private ArrayList<int[]> intTabNoGaps;
    private ArrayList<int[]> intTab;
    private ArrayList<int[]> intTabRC;
    private final String[] key;
    private final ArrayList<String[]> associativeList;
    ArrayList<char[]> alnTab;
    private static String Path, dirPath, SSZBINARY, RSCAPEBINARY;
    private double
            shannon = 0;
    private double SSZR_THRESHOLD = -3.0;        // alignments scoring below this will be kept (Z-score)
    private double percentAlignPower;
    private double bpCovary;
    private double totalBasePair;


    ScanItFast(ArrayList associativeList, String[] key, String Path, String dirPath,
               String SSZBINARY, String RSCAPEBINARY, boolean VERBOSE, boolean RSCAPE, boolean PRINTALL) {
        ScanItFast.Path = Path;
        ScanItFast.dirPath = dirPath;
        this.SSZBINARY = SSZBINARY;
        this.RSCAPEBINARY = RSCAPEBINARY;
        this.VERBOSE = VERBOSE;
        this.RSCAPE = RSCAPE;
        this.PRINTALL = PRINTALL;
        this.associativeList = associativeList;
        this.key = key;
    }

    public void run() {
        if (VERBOSE)
            System.out.println("- - -> Starting Scan");
        if (VERBOSE && alnTab.size() != associativeList.size()) {
            System.out.println(" #### Maf and AlnTab aren't same length");
            System.out.println(associativeList.size() + " " + alnTab.size());
        }

        int startPos = Integer.parseInt(key[1]);
//converts nucleotides to numbers

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
        Iterator iter = associativeList.iterator();

        intTabNoGaps = new ArrayList<>();
        intTab = new ArrayList<>();
        intTabRC = new ArrayList<>();
        ArrayList UniqueNames = new ArrayList();
        while (iter.hasNext()) {
            String[] line = (String[]) iter.next();
            String lineNoGap = line[1].replaceAll("[^ATCGUatcgu]", "").toUpperCase();
            String sequence = line[1].toUpperCase();
            int[] seqToInt = new int[sequence.length()];
            int[] seqToIntRC = new int[sequence.length()];
            int[] seqToIntNoGap = new int[lineNoGap.length()];
            for (int i = 0; i < sequence.length(); i++) {
                seqToInt[i] = letterMap.get(sequence.charAt(i));
                seqToIntRC[sequence.length() - i - 1] = letterMapRC.get(sequence.charAt(i));

            }
            for (int j = 0; j < lineNoGap.length(); j++) {
                int nuc = letterMap.get(lineNoGap.charAt(j));
                if (nuc != 4 || nuc != 5) {
                    seqToIntNoGap[j] = nuc;
                }
            }
            if (seqToIntNoGap.length > 0){
                intTabNoGaps.add(seqToIntNoGap);
                UniqueNames.add(line[0]);
                intTabRC.add(seqToIntRC);
                intTab.add(seqToInt);
            }
        }

        String[] nameTab = new String[UniqueNames.size()];
        nameTab = (String[]) UniqueNames.toArray(nameTab);
// first check for > 2 seqs
        int goodSeqs = intTab.size();

        if (goodSeqs <= 3) {
            if (VERBOSE)
                System.out.println("-> Not Enough seqs ");
            return;
        }
//check whether theres homo sapiens
        if (!(String.valueOf(UniqueNames.get(0)).startsWith("homo_sapiens"))){
            return;
        }

        // remove gappy sequences
        if (VERBOSE)
            System.out.println("- -> Gappy sequences");
        boolean[] keepMe = new boolean[intTab.size()];
        for (int seq = 0; seq != intTab.size(); seq++) {
            if (intTabNoGaps.get(seq).length >= intTab.get(seq).length * ((double) 50 / 100)) {
                keepMe[seq] = true;
            } else {
                keepMe[seq] = false;
                goodSeqs--;
            }

        }
        // System.out.println("this is the end of block");
        if (goodSeqs <= 3) {
            if (VERBOSE)
                System.out.println("-> Not Enough seqs in this window!");
            return;
        }
        // exit when human is shit
        if (!keepMe[0])
            return;

        // check for gap-only columns
        if (VERBOSE)
            System.out.println("- -> Gap only columns");
        boolean[] hasChars = new boolean[intTab.get(0).length];


/*********************************************************************
        					calculate stats						*
        *********************************************************************/
        if (VERBOSE)
            System.out.println("- - -> calculating statistics");
        double uniqueSeqs = goodSeqs;
        int outCols = intTab.get(0).length; //change last variable if CLUSTAL properties changes
        double[] stats = new double[6];
        double[] chars;
        double[] totalChars;
        boolean[] isNotUnique = new boolean[goodSeqs];

        double [] column = new double[(int) outCols];
        // calculate mpi
        for (int k = 0; k < outCols; k++) {
            double identicalNuc =0.0;
            double totalNuc =0.0;
            double [][] stats1 = new double[6][6];
            for (int i = 0; i != goodSeqs; i++) {
                for (int j = i + 1; j != goodSeqs; j++) {
                    stats1[intTab.get(i)[k]][intTab.get(j)[k]] += 1.0;

                }
            }
            for (int i=0; i<stats1.length; i++) {
                for (int j = 0; j < stats1[i].length; j++) {
                    totalNuc += stats1[i][j];
                    if (i == j) {
                        identicalNuc += stats1[i][j];
                    }
                }
            }
            column[k]= identicalNuc/totalNuc;
            double totalGaps = 0.0;
            for (int i=0; i< stats1.length; i++) {
                totalGaps += stats1[i][4] + stats1[i][5];
            }
            if (totalNuc != totalGaps){
                hasChars[k] = true;
            }

        }
        double sum=0.0;
        for(int v =0; v <column.length; v++){
            sum+=column[v];
        }
        double newMPI= sum/(column.length);




        totalChars = new double[]{0.0,0.0,0.0,0.0,0.0};
        //Convert sequence into numbers, 0 is a, 1 is b, 2 is c, 3 is g, 4 is N
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
            if ( keepMe[ seq ] ) {
                outAln[iterate] = nameTab[seq].substring(0, Math.min(nameTab[seq].length(), 20));
                outAlnRC[iterate] = nameTab[seq].substring(0, Math.min(nameTab[seq].length(), 20));
                for (int i = 0; i != 25 - Math.min(nameTab[seq].length(), 20); i++) {
                    outAln[iterate] = outAln[iterate] + " ";
                    outAlnRC[iterate] = outAlnRC[iterate] + " ";
                }
                for (int i = 0; i != intTab.get(0).length; i++) {
                    if (hasChars[i]) {
                        outAln[iterate] = outAln[iterate] + reverseMap.get(intTab.get(seq)[i]);
                    }
                    if (hasChars[intTabRC.get(0).length - i - 1]) {
                        outAlnRC[iterate] = outAlnRC[iterate] + reverseMap.get(intTabRC.get(seq)[i]);
                    }
                }
                outAln[iterate] = outAln[iterate] + "\n";
                outAlnRC[iterate] = outAlnRC[iterate] + "\n";
                iterate++;

            }
        }

        stats[0] = newMPI * 100;                                                                        // Mean Pairwise ID
        double var = 0.0;
        for (int i = 0; i < column.length; i++){
                var = var + (double) Math.pow(column[i] - newMPI, 2);
        }
        double standard = Math.sqrt(var/(column.length -1));
        stats[1] = standard;                        // Variance
        stats[2] = -1 * shannon / ((double) outCols);       // Normalized Shanon entropy
        stats[3] = 100 * (totalChars[2] + totalChars[3]) / (totalChars[0] + totalChars[1] + totalChars[2] + totalChars[3]);       // GC content
        stats[4] = 100 * totalChars[4] / (outCols * goodSeqs);// GAP content



        if(RSCAPE) {
        /*********************************************************************
        					R-scape scan & parse						*
        //*********************************************************************
        ///***************** 	R-scape scan & parse		******************/
        String stkName =  key[1] + "_" + key[2] + ".stk";
        String pathSpecName= dirPath + "/"+ key[1] + "_" + key[2];
        // Write Sequences to Stockholm Format
        File theDir = new File (pathSpecName);
        if (!theDir.exists()) {
            theDir.mkdir();
        }
        File stkFile = new File(theDir + "/" + stkName);

        try {
            BufferedWriter WriteStockholm = new BufferedWriter(new FileWriter( stkFile ));
            WriteStockholm.write("# STOCKHOLM 1.0\n") ;

            for (int i =0; i < associativeList.size(); i++){
                WriteStockholm.write(associativeList.get(i)[0] + "\t"+ associativeList.get(i)[1] +"\n");
            }

            String gcRef = "#=GC RF "+ key[6] + "\n";

            String  gcSScon =  "#=GC SS_cons " +key[7] + "\n";

            WriteStockholm.write(gcRef);
            WriteStockholm.write(gcSScon);
            WriteStockholm.write("//" +"\n");


            WriteStockholm.close() ;
        } catch (IOException Fuck) {
            if (VERBOSE)
                System.err.println("couldn't write stockholm file");
            Fuck.printStackTrace();
            stkFile.delete() ;
            return;
        }


        try {
            rScape(stkFile, RSCAPEBINARY, pathSpecName);
        } catch (IOException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }


        String endFileName = ".fold.power";
        File f = new File(pathSpecName);
        File[] matchingFiles = f.listFiles(new FilenameFilter() {
            public boolean accept(File dir, String name) {
                return name.endsWith(endFileName);
            }
        });

        String totalBasePairs = "";
        String expectedCovary = "";
        String observedCovary = "";


        try {

            Scanner scanner = new Scanner(matchingFiles[0]);

            //now read the file line by line...
            int lineNum = 0;
            int counter = 0;
            while (scanner.hasNextLine()) {
                String line = scanner.nextLine();
                lineNum++;

                if (line.startsWith("# BPAIRS expected to covary")) {
                    expectedCovary = line;

                } else if (line.startsWith("# BPAIRS observed to covary")) {
                    observedCovary = line;

                } else if (counter == 5) {
                    totalBasePairs = line;
                    counter++;

                } else if (line.startsWith("#")) {
                    counter++;

                }
            }

        } catch (FileNotFoundException e) {
            e.printStackTrace();
            return;
        }
        String[] expCov = expectedCovary.split(" ");
        String[] totalBP = totalBasePairs.split(" ");
        String[] obsCov = observedCovary.split(" ");
        double covExp = Double.parseDouble(expCov[5]);
        double totalBasePair = Double.parseDouble(totalBP[totalBP.length - 1]);
        double observCov = Double.parseDouble(obsCov[obsCov.length - 1]);
        double percentAlignPower;
        double bpCovary;
        if (totalBasePair != 0.0) {
            percentAlignPower = covExp / totalBasePair * 100;
            bpCovary = observCov / totalBasePair * 100;
        } else {
            percentAlignPower = 0.0;
            bpCovary = 0.0;
        }


        //Delete file with R-scape info
        String[] entries = theDir.list();
        for (String s : entries) {
            File currentFile = new File(theDir.getPath(), s);
            currentFile.delete();
        }
        theDir.delete();
    }



        // save BED coords
        if (VERBOSE)
            System.out.println("- -> Calculating BED coords ");
        String BedFile = key[0] + "\t";
        BedFile = BedFile + key[1] + "\t" + key[2] + "\t";

        BedFile = BedFile + goodSeqs + "\t" + ((double) (int) (10 * stats[0]) / 10) + "\t"      // MPI
                + ((double) (int) (10 * stats[1]) / 10) + "\t"            // STDEV
                + ((double) (int) (100 * stats[2]) / 100) + "\t"                // SHANNON
                + ((double) (int) (10 * stats[3]) / 10) + "\t"                  //      GC
                + ((double) (int) (10 * stats[4]) / 10);                     // GAPS




        if (VERBOSE)
            System.out.println("Pre SISSIz bed file: \n" + " " + BedFile);
        int random = (int) ((double) 10000 * Math.random());
        File Aln = new File(Path + "/" + BedFile.replaceAll("\t", "_") + ".aln." + random),    //
                AlnRC = new File(Path + "/" + BedFile.replaceAll("\t", "_") + "rc.aln." + random);  //
        // v v v v v v v v    INCLUSION STATS     v v v v v v v v v v v v v
       if ( stats[4] <= 75 && stats[0] > 50){
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
            } catch (IOException Fuck) {
                if (VERBOSE)
                    System.err.println("Arrgh... Couldn't write clustal file!");
                Fuck.printStackTrace();
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
        } catch (IOException Fuck) {
            Fuck.printStackTrace();
            System.err.println("ScanSSZ failed with ");
            for (int y = 0; y != goodSeqs; y++) {
                if (!isNotUnique[y]) {
                    System.err.println(outAln[y]);
                }
            }
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
                if (PRINTALL) {
                    System.out.println(FinalBedFile.replaceAll("_", "\t"));
                }
            } else {
                //write bed and rename alignment
                if(RSCAPE) {
                    System.out.println(FinalBedFile.replaceAll("_", "\t") + "\t" + percentAlignPower + "\t" +
                            bpCovary + "\t" + totalBasePair);
                }else {
                    System.out.println(FinalBedFile.replaceAll("_", "\t"));

                }
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
        } catch (IOException Fuck) {
            Fuck.printStackTrace();
            System.err.println("ScanSSZ failed in RC with ");
            for (int y = 0; y != goodSeqs; y++) {
                if (!isNotUnique[y]) {
                    System.err.println(outAln[y]);
                }
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
                if (PRINTALL) {
                    System.out.println(FinalBedFileRC.replaceAll("_", "\t"));
                }
            } else {
                if(RSCAPE){
                    System.out.println(FinalBedFileRC.replaceAll("_", "\t")+ "\t" + percentAlignPower
                            + "\t" + bpCovary + "\t" + totalBasePair);
                } else {
                    //write bedRC and rename alignment
                    System.out.println(FinalBedFileRC.replaceAll("_", "\t"));
                }
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
    					Rscape scan & parse						*
    //*********************************************************************/


    private void rScape(File stockholmFile, String RSCAPEBINARY, String path)throws IOException,
            InterruptedException {

        String cmd = RSCAPEBINARY + " --fold" + " --nofigures "
                + stockholmFile;
        Process process = Runtime.getRuntime().exec(cmd, null, new File(path));


        BufferedReader reader =
                new BufferedReader(new InputStreamReader(process.getInputStream()));
        while ((reader.readLine()) != null) {
        }
        process.waitFor();
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
            String SissizOutTab[] = new String[12];
            String Output, Error = "";
            String Command = SSZBINARY;

            Command = Command + " -j " + Path + "/" + BedFile.replaceAll("\t", "_") + ".aln." + id; // RIBOSUM scoring

            try {
                long now = System.currentTimeMillis();
                long timeoutInMillis = 1000L * 300;                          // max 5 minutes
                long finish = now + timeoutInMillis;
                // launch initial SISSIz call
                Process Sissiz = Runtime.getRuntime().exec(Command);
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
                        if (Output.length() > 6 && Output.startsWith("sissiz")) {
                            if (VERBOSE)
                                System.out.println(Output);
                            SissizOutTab = Output.split("\\s");
                            SissizOutTab[1] = "r";
                        }
                    }
                    SissizOut.close();
                }
                // rerun SISSIz if output is dodgy
                try {
                    /*>>>>>>*/
                    if (Double.parseDouble(SissizOutTab[7]) == 0
                            || (Math.abs(Double.parseDouble(SissizOutTab[9])) < 0.5)){ // variance of simulated alignments
                        if (VERBOSE)
                            System.out.println("SISSIz gave dodgy output... retrying with RIBOSUM");
                        SissizOutTab = new String[12];
                        now = System.currentTimeMillis();
                        finish = now + timeoutInMillis;
                        Sissiz = Runtime.getRuntime().exec(SSZBINARY + " -j " + Path + "/" + BedFile.replaceAll("\t", "_") + ".aln." + id);
                        if (VERBOSE)
                            System.out.println("Running " + SSZBINARY + " -j " + Path + "/" + BedFile.replaceAll("\t", "_") + ".aln." + id);
                        SissizErr = new BufferedReader(new InputStreamReader(Sissiz.getErrorStream()));
                        while (isAlive(Sissiz)) {
                            Thread.sleep(100);  // do we need this?
                            if (System.currentTimeMillis() > finish) {
                                if (VERBOSE)
                                    System.out.println("SISSIz failed to run within time");
                                SissizErr.close();
                                Sissiz.destroy();
                                return null;
                            } else {
                                // avoid infinity loop
                                if (SissizErr.ready()) {
                                    Error = SissizErr.readLine();

                                }
                            }
                        }
                        SissizErr.close();
                        if (SissizOutTab[0] == null) {
                            BufferedReader SissizOut = new BufferedReader(new InputStreamReader(Sissiz.getInputStream()));
                            while ((Output = SissizOut.readLine()) != null) {
                                if (Output.length() > 6 && Output.startsWith("sissiz")) {
                                    if (VERBOSE)
                                        System.out.println(Output + "  ...DONE !!\n");
                                    SissizOutTab = Output.split("\\s");
                                    SissizOutTab[1] = "r";
                                }
                            }
                            SissizOut.close();
                        }
                    }

                } catch (Exception Fuck) {
                    System.err.println("  with file = " + BedFile.replaceAll("\t", "_") + ".aln");
                    if (SissizOutTab[7] == null)
                        System.err.println("Null output data");
                    Fuck.printStackTrace();
                }
            } catch (Exception err) {
                System.out.println(" Caught Error!\n ----> " + Command + "\n  counter--> " );
                System.err.println("!!!Caught Error!\n ----> " + Command + "\n  counter--> " );
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
}

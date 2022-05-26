/**
  ECSFinder.java
  - - reads a set of maf files, feeds them to RNALalifold calculates stats,
  - - scans with SISSIz, outputs bed coordinates of high-confidence predictions
  Created by Vanda Gaonac'h-Lovejoy on 26/05/2022.
  Copyright 2022 All rights reserved.
*/
import java.util.*; import java.util.concurrent.* ; import java.io.*;
import java.lang.*;

public class ECSFinder {

    static int
            GAPS = 50,
            NTHREDS = 4;
    static boolean
            FILTER_ID = true, //for removing identical sequences
            VERBOSE = false,
            PRINTALL = false,
            REALIGN = false,
            RSCAPE = false;
    static String
            FILENAME = "",
            OUT_PATH = "",
            dirProgram = "",
            SSZBINARY = "~/SISSIz-master/src/SISSIz",
            ALIFOLDBINARY = "~/ViennaRNA-2.4.16/bin/RNALalifold",
            RSCAPEBINARY = "/usr/bin/R-scape";


    static double SSZR = -3.0;


    public static synchronized void main(String[] Args) throws IOException, InterruptedException {

        // variables
        String[] mafTabTemp;
        String[] TempTab;
        String Temp = "";
        String[] nameAlifold;
        int blockAln;
        // usage info
        if (Args.length == 0) {
            System.out.println("\n\t\t\t  version 1 \n" +

                    "\t SCAN MULTIPLE ALIGNMENTS FOR CONSERVED RNA STRUCTURES\n\n" +
                    "Reads a set of maf files, calculates stats, scans with SISSIz or RNAz, outputs bed coordonates of high-confidence predictions\n\n" +
                    "*** Known issues ***\n" +
                    "Usage:     java -jar ECSFinder.jar [options] -o output/directory -i input.maf (last parameter must be -i)\n\n" +
                    "Output: 	Two types of results are produced:" +
                    "           (1)  the multiple sequence alignments associated to significant predictions \n" +
                    "                are saved to files in the folder specified with the \"-o\" option.\n" +
                    "                File names correspond to their genomic coordinates in a .bed-compatible format. Ex:\n" +
                    "                     output/directory/chrX_12345_12500_80:75:23:14:8:z_300_+.aln\n" +
                    "           (2)  The genomic coordinates (.bed format) of ECSs are also written to the SDOUT\n" +
                    "                (see additional options below)." +
                    "output: Number chromosome\t start loci\t end loci\t number of species\t MPI\t standard deviation\t Normalized shannon entrpy" +
                            "GC % content \t Gap %content \t SISSIZ Z-score\n"+
                    "           ***  N.B. the score field corresponds to the SISSIz Z-score x-100\n\n" +
                    "Options:\n" +
                    "  -all             write out bed entires for all sampled alignments to STDOUT\n" +
                    "  -c     int       number of CPUs for calculations (default 4)\n" +
                    "  -g     int       max gap percentage of sequences for 2D prediction (default 50)\n" +
                    "  -sszr  double    report SISSIz+RIBOSUM hits below this Z-score (default -3)\n" +
                    "  -v               verbose (messy but detailed) output\n" +
                    " -r                R-scape to be used" );
            System.exit(0);
        }
        // get binary paths
        Process GetBinary = Runtime.getRuntime().exec("which RNALalifold");
        BufferedReader ReadBin = new BufferedReader(new InputStreamReader(GetBinary.getInputStream()));
        if ((ALIFOLDBINARY = ReadBin.readLine()) == null) {
            System.out.println("Please install RNALalifold and link it to your $PATH");
            System.exit(0);
        }

        GetBinary = Runtime.getRuntime().exec("which SISSIz") ;
        ReadBin = new BufferedReader(new InputStreamReader(GetBinary.getInputStream()));
        if ( (SSZBINARY = ReadBin.readLine() ) == null ) {
            System.out.println("Please install SISSIz (version 2.0), and link it to your $PATH" );
            System.exit(0);
        }

        GetBinary = Runtime.getRuntime().exec("which R-scape") ;
        ReadBin = new BufferedReader(new InputStreamReader(GetBinary.getInputStream()));

        ReadBin.close();

        // parse arguments
        for (int i = 0; i != Args.length; i++) {
            if (Args[i].equals("-c")) {  // Threads
                NTHREDS = Integer.parseInt(Args[i + 1]);
                i++;
            } else if (Args[i].equals("-f")) {  // Threads
                FILTER_ID = false;
                i++;
            } else if (Args[i].equals("-g")) {  // gap content
                GAPS = Integer.parseInt(Args[i + 1]);
                i++;
            } else if (Args[i].equals("-o")) {//output directory
                OUT_PATH = System.getProperty("user.dir") + "/"+ Args[i + 1];
                if (!(new File(OUT_PATH)).isDirectory())
                    (new File(OUT_PATH)).mkdirs();
                if (VERBOSE)
                    System.out.println("writing alignments to directory " + OUT_PATH);
                dirProgram = System.getProperty("user.dir");
                i++;
            } else if (Args[i].equals("-all")) { // step size
                PRINTALL = true;
                i++;
            } else if (Args[i].equals("-v")) { // verbose output
                VERBOSE = true;
            }else if (Args[i].equals("-r")){ //Rscape analysis
                RSCAPE = true;
            } else if (Args[i].equals("-mafft")) { // realign
                REALIGN = true;
            } else if (Args[i].equals("-i")) {
                i++;
                FILENAME = Args[i].substring(Args[i].lastIndexOf("/") + 1);
                nameAlifold = FILENAME.split("\\.");
                //parse out individual alignment blocks from a multi maf file
                int lineCount = 0;
                BufferedReader ReadFile = new BufferedReader(new FileReader(Args[i]));
                String Line;
                while ((Line = ReadFile.readLine()) != null)   // count lines for array
                    if (Line.length() > 1 && Line.charAt(0) != '#')
                        lineCount++;

                if (VERBOSE)
                    System.out.println("Read " + (lineCount - 1) + " sequences from file " + FILENAME);
                ReadFile.close();

                /************************************************************************
                 ****   RNALalifold       ****
                 ************************************************************************/
                String cmd = ALIFOLDBINARY + " --id-prefix=alifold" + " --noLP" + " --maxBPspan=300" + " --ribosum_scoring"
                        + " --aln-stk " + Args[Args.length - 1];
                String Path2 =  OUT_PATH + "/stockholm" + nameAlifold[nameAlifold.length - 1];
                File stockholmFolder = new File(Path2);
                if (!stockholmFolder.exists()) {
                    stockholmFolder.mkdir();
                }
                 executeCommand(cmd, nameAlifold);

                ReadFile = new BufferedReader(new FileReader(Args[i]));
                // fill in array from file
                blockAln = 0;
                /************************************************************************
                 ****   This stuff is messy, but avoids problems at last block       ****
                 ************************************************************************/

                ExecutorService MultiThreads = Executors.newFixedThreadPool(NTHREDS);
                List<Future<Runnable>> futures;
                readBlocks:
                while ((Line = ReadFile.readLine()) != null ) {
                    if (Line.length() >= 1 && Line.charAt(0) != '#') {
                        if (Line.charAt(0) == 'a') {
                            blockAln++;
                            continue readBlocks;
                        } else if (Line.substring(0, 1).equals("s")) {
                            Temp = Temp + Line + "@";
                        }

                    } else if ((Temp.split("@").length <= 2) && Line.equals("")){
                        Temp = "";
                    }else if ((Temp.split("@").length >= 3) && Line.equals("")) { // at least 3 sequences
                        TempTab = Temp.split("@");
                        Temp = "";

                        ArrayList<String[]> associativeList = new ArrayList<>();
                        mafTabTemp = TempTab[0].split("\\s+");
                        futures = new ArrayList<Future<Runnable>>();
                        // add Path Flag ? < < < < < < < < < < < < < < < < < < < < < < < < < < < < <
                        String Path = OUT_PATH + "/aln/" + mafTabTemp[1].substring(mafTabTemp[1].lastIndexOf(".") + 1);
                        if (!(new File(Path)).isDirectory())
                            (new File(Path)).mkdirs();

                        int numbZeros = 0;
                        String gcReference = "";
                        String gcSScons = "";
                        String lastDigit = "";
                        String finalName = "";
                        lastDigit += String.valueOf(blockAln);
                        numbZeros = 4 - lastDigit.length();
                        int x = 0;
                        while (x < numbZeros) {
                            finalName += String.valueOf(0);
                            x++;
                        }
                        finalName += lastDigit;

                        File file = new File( OUT_PATH + "/stockholm" +
                                nameAlifold[nameAlifold.length - 1] + "/alifold_"
                                + finalName + ".stk");

                        if (file.exists()){
                            try {

                                System.out.println("Stockholm file: " +
                                        file.getAbsolutePath().substring(file.getAbsolutePath().lastIndexOf("/")
                                                + 1));
                                BufferedReader reader = new BufferedReader(new FileReader(file));
                                String currentLine = reader.readLine();

                                String[] arrayName = new String[5];
                                int lengthAln = 0;

                                readAlns:

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
                                        continue readAlns;

                                    } else if ( !(currentLine.startsWith("//"))) {

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

                                        if (Integer.parseInt(arrayLociChrm[2])-Integer.parseInt(arrayLociChrm[1]) < 50){
                                            currentLine = reader.readLine();
                                            continue readAlns;
                                        }
                                        ScanItFast aln = new ScanItFast(associativeList,
                                                arrayLociChrm, Path, OUT_PATH,
                                                SSZBINARY, RSCAPEBINARY, VERBOSE, RSCAPE, PRINTALL);
                                        aln.setSszR(SSZR);

                                        Future f = MultiThreads.submit(aln);
                                        futures.add(f);
                                    }
                                    currentLine = reader.readLine();
                                    continue readAlns;
                                }


                                reader.close();


                             //     file.delete();


                            } catch (FileNotFoundException e) {
                                e.printStackTrace();
                            }

                    } else {
                            System.out.println("File: "+ file.getAbsolutePath().substring(file.getAbsolutePath().lastIndexOf("/")
                                    + 1) + " does not exist");
                        }
                                for (Future<Runnable> g : futures) {
                                    try {
                                        g.get();


                                    } catch (Exception Fuck) {
                                        System.err.println("MultiThreads took too long!  OOps!");
                                        Fuck.printStackTrace();
                                    }
                                }

                            }

                    }
                //Delete file with stockholm info
                String[] entries = stockholmFolder.list();
                if (entries.length > 0) {
                    for (String s : entries) {
                        File currentFile = new File(stockholmFolder.getPath(), s);
                        currentFile.delete();
                    }
                }
                stockholmFolder.delete();

                ReadFile.close();
                MultiThreads.shutdown();
                MultiThreads.awaitTermination(60 * 10L, TimeUnit.SECONDS);

            }
        }




    }

    /*********************************************************************
    					Alignment generator				*
    //*********************************************************************/



    private static int[] getRealCoordinates (int start, String[] mafCord, String motifHuman){
        int [] cordFinal;
        int [] cordFinalPlus1= new int[2];
        String withoutGap= mafCord[6].substring(0, start);

        int nuc = withoutGap.replaceAll("-", "").length();

        if (mafCord[4].equals("-")){
            int lociEnd = (Integer.parseInt(mafCord[5] ) + 1  - (Integer.parseInt(mafCord[2]) + nuc)) + 1 ;

            int lociStart = (lociEnd - motifHuman.replaceAll("-", "").length());

            cordFinal = new int[]{lociStart, lociEnd};
        } else {
            int lociStart = (Integer.parseInt(mafCord[2]) + nuc) ;
            int lociEnd = lociStart + (motifHuman.replaceAll("-", "").length());
            cordFinal = new int[]{lociStart , lociEnd};
        }
        cordFinalPlus1[0]= cordFinal[0];
        cordFinalPlus1[1]= cordFinal[1] -1;

        return cordFinalPlus1;
    }


    private static void executeCommand(final String command, String[] nameAlifold) {
        String Path = OUT_PATH + "/stockholm" + nameAlifold[nameAlifold.length - 1];
        System.out.println("Executing command " + command);
        try {

            Process process = Runtime.getRuntime().exec(command, null, new File(Path));
            InputStream error = process.getInputStream();
            for(int i=0; i<error.available(); i++){
                System.out.println(""+ error.read());
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

}


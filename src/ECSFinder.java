/*
  ECSFinder.java
  - - reads a set of maf files, feeds them to RNALalifold calculates stats,
  - - scans with SISSIz, outputs bed coordinates of high-confidence predictions
  Created by Vanda Gaonac'h-Lovejoy on 26/05/2022.
  Copyright 2022 All rights reserved.
*/

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

public class ECSFinder {

    static int
            GAPS = 50,
            NTHREDS = 4;
    static boolean
            VERBOSE = false;
    static String
            FILENAME = "",
            OUT_PATH = "",
            dirProgram = "",
            SSZBINARY = "~/SISSIz-master/src/SISSIz",
            ALIFOLDBINARY = "~/ViennaRNA-2.4.16/bin/RNALalifold";


    static double SSZR = -3.0;


    public static synchronized void main(String[] Args) throws IOException, InterruptedException {

        // variables
        String[] mafTabTemp;
        String[] TempTab;
        StringBuilder Temp = new StringBuilder();
        String[] nameAlifold;
        int blockAln;
        // usage info
        if (Args.length == 0) {
            System.out.println("\n\t\t\t  version 1.0 \n" +
                    " ________    ______   ______    ________  _                 __                \n"+
                    "|_   __  | .' ___  |.' ____ \\  |_   __  |(_)               |  ]               \n"+
                    "  | |_ \\_|/ .'   \\_|| (___ \\_|   | |_ \\_|__   _ .--.   .--.| | .---.  _ .--.  \n"+
                    "  |  _| _ | |        _.____`.    |  _|  [  | [ `.-. |/ /'`\\' |/ /__\\\\[ `/'`\\] \n"+
                    " _| |__/ |\\ `.___.'\\| \\____) |  _| |_    | |  | | | || \\__/  || \\__., | |     \n"+
                    "|________| `.____ .' \\______.' |_____|  [___][___||__]'.__.;__]'.__.'[___]    \n"+




            "\t SCAN MULTIPLE ALIGNMENTS FOR CONSERVED RNA STRUCTURES\n\n" +
                    "Reads a set of maf files, calculates stats, scans with SISSIz , outputs bed coordonates of high-confidence predictions\n\n" +
                    "*** Known issues ***\n" +
                    "Usage:     java  ECSFinder [options] -o output/directory -i input.maf (last parameter must be -i)\n\n" +
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
                    "  -c     int       number of CPUs for calculations (default 4)\n" +
                    "  -g     int       max gap percentage of sequences for 2D prediction (default 50)\n" +
                    "  -sszr  double    report SISSIz+RIBOSUM hits below this Z-score (default -3.0)\n" +
                    "  -v               verbose (messy but detailed) output\n");
            System.exit(0);
        }
        // get binary paths
        ProcessBuilder pbCmd = new ProcessBuilder("which","RNALalifold");
        Process pbCmdRNAL = pbCmd.start();
        BufferedReader ReadBin = new BufferedReader(new InputStreamReader(pbCmdRNAL.getInputStream()));
        if ((ALIFOLDBINARY = ReadBin.readLine()) == null) {
            System.out.println("Please install RNALalifold and link it to your $PATH");
            System.exit(0);
        }

        ProcessBuilder pbCmdSISSIz = new ProcessBuilder("which", "SISSIz") ;
        Process pbCmdSISSIz2 = pbCmdSISSIz.start();
        ReadBin = new BufferedReader(new InputStreamReader(pbCmdSISSIz2.getInputStream()));
        if ( (SSZBINARY = ReadBin.readLine() ) == null ) {
            System.out.println("Please install SISSIz (version 2.0), and link it to your $PATH" );
            System.exit(0);
        }


        ReadBin.close();

        // parse arguments
        for (int i = 0; i != Args.length; i++) {
            switch (Args[i]) {
                case "-c":   // Threads
                    NTHREDS = Integer.parseInt(Args[i + 1]);
                    i++;
                    break;
                case "-g":   // gap content
                    GAPS = Integer.parseInt(Args[i + 1]);
                    i++;
                    break;
                case "-o": //output directory
                    OUT_PATH = System.getProperty("user.dir") + "/" + Args[i + 1];
                    if (!(new File(OUT_PATH)).isDirectory())
                        (new File(OUT_PATH)).mkdirs();
                    if (VERBOSE)
                        System.out.println("writing alignments to directory " + OUT_PATH);
                    dirProgram = System.getProperty("user.dir");
                    i++;

                    break;
                case "-v":  // verbose output
                    VERBOSE = true;
                    break;
                case "-sszr":  // verbose output
                    SSZR=Double.valueOf(Args[i + 1]);
                    break;
                case "-i":
                    i++;
                    FILENAME = Args[i].substring(Args[i].lastIndexOf("/") + 1);

                    if (Args[i].split("/").length==1){
                        System.out.println("Please enter the full path of the input file");
                        System.exit(0);
                    }
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


                    String Path2 = OUT_PATH + "/stockholm" + nameAlifold[nameAlifold.length - 1];
                    File stockholmFolder = new File(Path2);
                    if (!stockholmFolder.exists()) {
                        stockholmFolder.mkdir();
                    }
                    //executeCommand(Args[Args.length - 1],nameAlifold);

                    ReadFile = new BufferedReader(new FileReader(Args[i]));
                    // fill in array from file
                    blockAln = 0;
                    /************************************************************************
                     ****   This stuff is messy, but avoids problems at last block       ****
                     ************************************************************************/

                    ExecutorService MultiThreads = Executors.newFixedThreadPool(NTHREDS);
                    List<Future<?>> futures;
                    while ((Line = ReadFile.readLine()) != null) {
                        if (Line.length() >= 1 && Line.charAt(0) != '#') {
                            if (Line.charAt(0) == 'a') {
                                blockAln++;
                            } else if (Line.charAt(0) == 's') {
                                Temp.append(Line).append("@");
                            }

                        } else if ((Temp.toString().split("@").length <= 2) && Line.equals("")) {
                            Temp = new StringBuilder();
                        } else if ((Temp.toString().split("@").length >= 3) && Line.equals("")) { // at least 3 sequences
                            TempTab = Temp.toString().split("@");
                            Temp = new StringBuilder();

                            ArrayList<String[]> associativeList = new ArrayList<>();
                            mafTabTemp = TempTab[0].split("\\s+");
                            futures = new ArrayList<>();
                            String Path = OUT_PATH + "/aln/" + mafTabTemp[1].substring(mafTabTemp[1].lastIndexOf(".") + 1);
                            if (!(new File(Path)).isDirectory())
                                (new File(Path)).mkdirs();

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
                                            ScanItFast aln = new ScanItFast(associativeList, arrayLociChrm, Path,
                                                    SSZBINARY, VERBOSE);
                                            aln.setSszR(SSZR);
                                            aln.setGap(GAPS);

                                            Future<?> f =  MultiThreads.submit(aln);
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
                    /*//Delete file with stockholm info
                    String[] entries = stockholmFolder.list();
                    assert entries != null;
                    if (entries.length > 0) {
                        for (String s : entries) {
                            File currentFile = new File(stockholmFolder.getPath(), s);
                            currentFile.delete();
                        }
                    }
                    stockholmFolder.delete();
*/
                    ReadFile.close();
                    MultiThreads.shutdown();
                    MultiThreads.awaitTermination(60 * 10L, TimeUnit.SECONDS);

                    break;
            }
        }



    }

    /*********************************************************************
    					Get coordinates				*
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


    private static void executeCommand(final String file, String[] nameAlifold) {
        String Path = OUT_PATH + "/stockholm" + nameAlifold[nameAlifold.length - 1];

        try {

            ProcessBuilder pb = new ProcessBuilder(ALIFOLDBINARY, "--id-prefix=alifold", "--noLP",
                    "--maxBPspan=300","--ribosum_scoring", "--aln-stk", file );
            String cmd =ALIFOLDBINARY + " --id-prefix=alifold"+ " --noLP"+
                    " --maxBPspan=300"+" --ribosum_scoring"+ " --aln-stk "+ file;
            if (VERBOSE)
            System.out.println("Executing command " + cmd);

            pb.directory(new File(Path));
            Process process = pb.start();
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


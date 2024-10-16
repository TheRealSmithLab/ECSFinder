
import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class ECSFinder {

    static int GAPS = 50, NTHREDS = 4;
    static boolean VERBOSE = false, MAFFT =false;
    static String FILENAME = "", OUT_PATH = "", dirProgram = "";
    static String SSZBINARY = "", ALIFOLDBINARY = "",
            RNAALIFOLD="",R= "", RSCAPE="", MAFFTBINARY="";
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
                    " -mafft realign aln using MAFFT (default FALSE)\n"+
                    " -v verbose (messy but detailed) output\n");
            System.exit(0);
        }



        // parse arguments
        parseArguments(args);
        // get binary paths
        setBinaryPaths();
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
        if(MAFFT){
            MAFFTBINARY = getBinaryPath("mafft-ginsi");
        }
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
                    NTHREDS = Integer.parseInt(args[++i]);
                    break;
                case "-g":
                    GAPS = Integer.parseInt(args[++i]);
                    break;
                case "-o":
                    OUT_PATH = System.getProperty("user.dir") + "/" + args[++i];
                    createDirectory(OUT_PATH);
                    dirProgram = System.getProperty("user.dir");
                    break;
                case "-v":
                    VERBOSE = true;
                    break;
                case "-mafft":
                    MAFFT = true;
                    break;
                case "-sszr":
                    SSZR = Double.parseDouble(args[++i]);
                    break;
                case "-i":
                    FILENAME = args[++i];
                    break;
                default:
                    System.err.println("Invalid argument: " + args[i]);
                    printUsageAndExit();
            }
        }

        // Additional argument validation
        if (FILENAME.isEmpty()) {
            System.err.println("Input MAF file (-i) is required.");
            printUsageAndExit();
        }
        if (OUT_PATH.isEmpty()) {
            System.err.println("Output directory (-o) is required.");
            printUsageAndExit();
        }
    }

    private static void printUsageAndExit() {
        System.out.println("Usage: java ECSFinder [options] -o output/directory -i input.maf");
        System.exit(1);
    }


    private static void createDirectory(String path) {
        File dir = new File(path);
        if (!dir.isDirectory()) {
            dir.mkdirs();
        }
    }

    private static void preprocessMafFiles() {
        File inputDir = new File(FILENAME);
        if (!inputDir.isDirectory()) {
            System.err.println("Error: Provided path is not a directory or does not exist.");
            System.exit(1);
        }

        File[] mafFiles = inputDir.listFiles((dir, name) -> name.endsWith(".maf"));
        if (mafFiles == null || mafFiles.length == 0) {
            System.err.println("Error: No .maf files found in the directory: " + FILENAME);
            System.exit(1);
        }

        String[] mafFilePaths = Arrays.stream(mafFiles).map(File::getAbsolutePath).toArray(String[]::new);

        try {
            MergeNFilter mergeNFilter = new MergeNFilter();
            mergeNFilter.process(mafFilePaths, OUT_PATH);
        } catch (IOException e) {
            System.err.println("Error processing MAF files: " + e.getMessage());
            e.printStackTrace();
            System.exit(1);
        }
    }


    private static void runRNALalifoldAndProcessResults() throws IOException, InterruptedException {
        // Extract file path
        String inputFile = OUT_PATH + "/output.maf";
        String stockholmFolderPath = OUT_PATH + "/stockholm";
        File stockholmFolder = createStockholmFolder(stockholmFolderPath);

        // Handle MAFFT or direct MAF file processing
        if (!MAFFT) {
            runRNALalifold(inputFile);
        } else {
            processWithMafft(inputFile);  // MAFFT handling
        }

        // Process MAF file for alignment blocks
        try {
            processAlignmentBlocks(inputFile);
        } catch (ExecutionException e) {
            System.err.println("Error in executing task: " + e.getMessage());
            e.printStackTrace();
        }

        // Clean up the Stockholm folder
        cleanUpFolder(stockholmFolder);
    }


    private static File createStockholmFolder(String stockholmFolderPath) {
        File stockholmFolder = new File(stockholmFolderPath);
        if (!stockholmFolder.exists()) {
            stockholmFolder.mkdir();
        }
        return stockholmFolder;
    }

    private static void processWithMafft(String inputFile) throws IOException, InterruptedException {
        // Convert MAF to FASTA and realign with MAFFT
        convertMafToSeparateFastas(inputFile, OUT_PATH);
        realignFastaFilesWithMafft();  // Realign all FASTA files in the specified directory
        // Process all realigned FASTA files
        File outputDir = new File(OUT_PATH + "/output_fasta_dir");
        File[] fastaFilesReAln = outputDir.listFiles((dir, name) -> name.endsWith("_realigned.fasta"));

        if (fastaFilesReAln != null) {
            // Iterate through each realigned FASTA file
            for (int i = 0; i < fastaFilesReAln.length; i++) {
                String realignedFilePath = fastaFilesReAln[i].getAbsolutePath();
                // Pass the block number 'i + 1' as the second argument
                runRNALalifold(realignedFilePath);  // Assign a block number for each FASTA file
            }
        }
        cleanUpFolder(outputDir);
    }


    private static void processAlignmentBlocks(String inputFile) throws IOException, InterruptedException, ExecutionException {
        BufferedReader readFile = new BufferedReader(new FileReader(inputFile));
        int blockAln = 0;
        StringBuilder temp = new StringBuilder();
        ExecutorService multiThreads = Executors.newFixedThreadPool(NTHREDS);
        List<Future<?>> futures = new ArrayList<>();
        path = new File(OUT_PATH + "/aln/");
        if (!path.isDirectory()) {
            path.mkdirs();  // Create the directory if it doesn't exist
        }

        String line;
        while ((line = readFile.readLine()) != null) {
            if (line.length() >= 1 && line.charAt(0) != '#') {
                if (line.charAt(0) == 'a') {
                    blockAln++;
                } else if (line.charAt(0) == 's') {
                    temp.append(line).append("@");
                }
            } else if ((temp.toString().split("@").length <= 2) && line.equals("")) {
                // If the block has fewer than 3 sequences, print a warning and skip the block
                System.out.println("Warning: Alignment block " + blockAln + " has fewer than 3 sequences. Skipping this block.");
                temp = new StringBuilder();
            } else if ((temp.toString().split("@").length >= 3) && line.equals("")) { // At least 3 sequences
                processAlignmentBlock(temp.toString(), blockAln,  futures, multiThreads);
                temp = new StringBuilder();  // Reset for the next block
            }
        }

        // Wait for all threads to complete and handle exceptions
        try {
            for (Future<?> future : futures) {
                future.get();  // This line can throw ExecutionException or InterruptedException
            }
        } catch (ExecutionException e) {
            System.err.println("Error in task execution: " + e.getMessage());
            throw e;  // Rethrow or handle as needed
        } catch (InterruptedException e) {
            System.err.println("Task was interrupted: " + e.getMessage());
            throw e;  // Rethrow or handle as needed
        }

        multiThreads.shutdown();
        multiThreads.awaitTermination(60 * 10L, TimeUnit.SECONDS);

        readFile.close();
    }

    private static void processAlignmentBlock(String block, int blockAln, List<Future<?>> futures, ExecutorService multiThreads) throws IOException {
        ArrayList<String[]> associativeList = new ArrayList<>();
        String[] mafTabTemp = block.split("@")[0].split("\\s+");

        // Generate Stockholm file names based on block number
        String finalName = getBlockName(blockAln);
        File stockholmFile = new File(OUT_PATH + "/stockholm"  + "/alifold_" + finalName + ".stk");

        if (stockholmFile.exists()) {
            processStockholmFile(stockholmFile, mafTabTemp, associativeList, futures, multiThreads);
        } else {
            System.out.println("File: " + stockholmFile.getName() + " does not exist");
        }
    }

    private static void processStockholmFile(File stockholmFile, String[] mafTabTemp, ArrayList<String[]> associativeList, List<Future<?>> futures, ExecutorService multiThreads) {
        try (BufferedReader reader = new BufferedReader(new FileReader(stockholmFile))) {
            String currentLine;
            String[] arrayName = new String[5];
            String gcReference = "", gcSScons = "";

            while ((currentLine = reader.readLine()) != null) {
                if (currentLine.startsWith("#=GF ID ")) {
                    arrayName = currentLine.split("[_.]");
                    associativeList = new ArrayList<>();
                } else if (currentLine.startsWith("#=GC RF")) {
                    gcReference = extractValue(currentLine);
                } else if (currentLine.startsWith("#=GC SS_cons")) {
                    gcSScons = extractValue(currentLine);
                } else if (!currentLine.startsWith("#") && !currentLine.equals("") && !currentLine.startsWith("//")) {
                    associativeList.add(processSpeciesLine(currentLine));
                }

                if ((!associativeList.isEmpty()) && currentLine.startsWith("//")) {
                    processMotif(mafTabTemp, arrayName, associativeList, gcReference, gcSScons, futures, multiThreads);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void processMotif(String[] mafTabTemp, String[] arrayName, ArrayList<String[]> associativeList, String gcReference, String gcSScons, List<Future<?>> futures, ExecutorService multiThreads) {
        int[] cordMotif = getRealCoordinates(Integer.parseInt(arrayName[3]), mafTabTemp, associativeList.get(0)[1]);
        String loci = Arrays.toString(cordMotif);
        String chrom = mafTabTemp[1].substring((mafTabTemp[1].lastIndexOf(".") + 1));
        String lociChrm = chrom + ", " + loci.substring(1, loci.length() - 1) + ", " + mafTabTemp[4] + ", " + arrayName[3] + ", " + arrayName[4] + ", " + gcReference + ", " + gcSScons;
        String[] arrayLociChrm = lociChrm.split(", ");

        if (Integer.parseInt(arrayLociChrm[2]) - Integer.parseInt(arrayLociChrm[1]) < 50) {
            return;  // Skip short alignments
        }

        // Create the ScanItFast task and submit it to the thread pool
        ScanItFast aln = new ScanItFast(associativeList, arrayLociChrm, path, SSZBINARY, VERBOSE);
        aln.setSszR(SSZR);
        aln.setGap(GAPS);
        Future<?> future = multiThreads.submit(aln);
        futures.add(future);
    }

    private static String getBlockName(int blockAln) {
        String lastDigit = String.valueOf(blockAln);
        int numbZeros = 4 - lastDigit.length();
        StringBuilder finalName = new StringBuilder();
        for (int i = 0; i < numbZeros; i++) {
            finalName.append('0');
        }
        finalName.append(lastDigit);
        return finalName.toString();
    }

    private static String extractValue(String line) {
        String[] lineReference = line.split(" ");
        return lineReference[lineReference.length - 1];
    }

    private static String[] processSpeciesLine(String line) {
        String[] species = line.split(" ", 2);
        species[1] = species[1].trim();
        return species;
    }


    private static void cleanUpFolder(File stockholmFolder) {
        String[] entries = stockholmFolder.list();
        if (entries != null && entries.length > 0) {
            for (String entry : entries) {
                File currentFile = new File(stockholmFolder.getPath(), entry);
                currentFile.delete();
            }
        }
        stockholmFolder.delete();
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


    /**
     * Converts a MAF file to separate FASTA files, one per alignment block.
     *
     * @param mafFilePath   Path to the input MAF file.
     * @param outputDirPath Directory to output separate FASTA files.
     */
    public static void convertMafToSeparateFastas(String mafFilePath, String outputDirPath) {
        String fastaOutput = "output_fasta_dir";  // Directory where the FASTA files will be saved
        int blockCount = 0;  // To track the block number
        BufferedWriter writer = null;

        try (BufferedReader reader = new BufferedReader(new FileReader(mafFilePath))) {
            // Create the directory if it does not exist
            Files.createDirectories(Paths.get(outputDirPath + "/" + fastaOutput));

            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();  // Trim whitespace

                if (line.startsWith("#") || line.isEmpty()) {
                    // Skip comments and empty lines
                    continue;
                }

                if (line.startsWith("a")) {
                    // A new alignment block starts
                    blockCount++;

                    // Close previous block's writer (if any)
                    if (writer != null) {
                        writer.close();
                    }

                    // Create a new FASTA file for this block
                    String blockFileName = outputDirPath + "/" + fastaOutput + "/block_" + blockCount + ".fasta";
                    writer = new BufferedWriter(new FileWriter(blockFileName));
                } else if (line.startsWith("s")) {
                    // Sequence line; extract relevant info
                    String[] tokens = line.split("\\s+");
                    String sequenceId = tokens[1];  // Sequence ID
                    String sequence = tokens[tokens.length - 1];  // Aligned sequence with gaps

                    // Write the FASTA formatted sequence to the current block's FASTA file
                    writer.write(">" + sequenceId + "\n");
                    writer.write(sequence + "\n");
                }
            }

            // Close the last block's writer
            if (writer != null) {
                writer.close();
            }

        } catch (IOException e) {
            System.err.println("Error processing the MAF file: " + e.getMessage());
        }
    }


    /**
     * Realigns sequences using MAFFT and writes the output directly to a FASTA file.
     * @param inputFilePath    The path to the input FASTA file.
     * @param outputFilePath   The path to the output FASTA file.
     * @throws IOException           If an I/O error occurs.
     * @throws InterruptedException  If the MAFFT process is interrupted.
     */
    public static void realignSequences(File inputFilePath, File outputFilePath) throws IOException, InterruptedException {
        String realignedFilePath = inputFilePath.getAbsolutePath().replace(".fasta", "_realigned.fasta");
        // Strip gaps from the input FASTA sequences before realignment
        File tempInputFile = new File(outputFilePath+"/temp_gap_stripped.fasta");

        try (BufferedReader reader = new BufferedReader(new FileReader(inputFilePath));
             BufferedWriter tempWriter = new BufferedWriter(new FileWriter(tempInputFile))) {

            String line;
            StringBuilder sequence = new StringBuilder();
            String header = "";

            while ((line = reader.readLine()) != null) {
                if (line.startsWith(">")) {
                    // Write the previous sequence if any (after stripping gaps)
                    if (sequence.length() > 0) {
                        tempWriter.write(sequence.toString().replace("-", "") + "\n");
                        sequence.setLength(0);  // Clear the sequence for the next one
                    }
                    // Write the header to the temp file
                    header = line;
                    tempWriter.write(header + "\n");
                } else {
                    // Accumulate sequence lines (gaps will be removed later)
                    sequence.append(line.trim());
                }
            }
            // Write the last sequence
            if (sequence.length() > 0) {
                tempWriter.write(sequence.toString().replace("-", "") + "\n");
            }
        }

        //Run MAFFT on the gap-stripped sequences
        List<String> command = Arrays.asList(MAFFTBINARY, "--quiet", "--thread", String.valueOf(NTHREDS), tempInputFile.getAbsolutePath());

        ProcessBuilder pb = new ProcessBuilder(command);
        Process mafftProcess = pb.start();

        // Write the realigned sequences to the final output FASTA file
        try (BufferedReader reAligned = new BufferedReader(new InputStreamReader(mafftProcess.getInputStream()));
             BufferedWriter writer = new BufferedWriter(new FileWriter(realignedFilePath))) {

            StringBuilder sequence = new StringBuilder();
            String speciesReal = "";
            String line;

            while ((line = reAligned.readLine()) != null) {
                if (line.startsWith(">")) {
                    // If we encounter a new sequence header, write the previous sequence to the file
                    if (!speciesReal.isEmpty()) {
                        writer.write(">" + speciesReal + "\n");
                        writer.write(sequence.toString().toUpperCase() + "\n");
                        sequence = new StringBuilder();  // Reset for the next sequence
                    }
                    speciesReal = line.substring(1);  // Remove '>' from the header line
                } else {
                    // Append sequence lines
                    sequence.append(line);
                }
            }

            // Write the last sequence to the FASTA file
            if (!speciesReal.isEmpty()) {
                writer.write(">" + speciesReal + "\n");
                writer.write(sequence.toString().toUpperCase() + "\n");
            }
        }

        // Wait for the MAFFT process to finish and check the exit code
        int exitCode = mafftProcess.waitFor();
        if (exitCode != 0) {
            throw new IOException("MAFFT exited with error code: " + exitCode);
        }

        //  Clean up temporary files (optional)
        if (tempInputFile.exists()) {
            tempInputFile.delete();
        }
    }


    /**
     * Saves the realigned sequences to a new FASTA file.
     *
     * @param realignedSequences The map of realigned sequences to save.
     * @param outputFilePath     The path to the output FASTA file.
     * @throws IOException
     */
    public static void saveRealignedSequences(Map<String, String> realignedSequences, String outputFilePath) throws IOException {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFilePath))) {
            for (Map.Entry<String, String> entry : realignedSequences.entrySet()) {
                writer.write(">" + entry.getKey() + "\n");
                writer.write(entry.getValue() + "\n");
            }
            System.out.println("Realigned sequences saved to: " + outputFilePath);
        }
    }
    /**
     * Realigns all FASTA files in the given directory using MAFFT.
     *
     * @throws IOException           If an I/O error occurs.
     * @throws InterruptedException  If the MAFFT process is interrupted.
     */
    private static void realignFastaFilesWithMafft() throws IOException, InterruptedException {
        File outputDir = new File(OUT_PATH + "/output_fasta_dir");
        File[] fastaFiles = outputDir.listFiles((dir, name) -> name.endsWith(".fasta"));

        if (fastaFiles != null) {
            for (File fastaFile : fastaFiles) {

                realignSequences(fastaFile, outputDir);
            }
        }
    }
    /**
     * Executes RNALalifold on the given input file (MAF or FASTA) and writes the result to a unique Stockholm file.
     * For MAF, it automatically generates files like alifold_0001.stk, alifold_0002.stk, etc.
     * For FASTA, we manually handle unique filenames by appending block numbers.
     *
     * @param inputFilePath  Path to the input alignment file (either MAF or FASTA)
     * @throws IOException           If an I/O error occurs.
     * @throws InterruptedException  If the process is interrupted.
     */
    private static void runRNALalifold(String inputFilePath) throws IOException, InterruptedException {

            String outputFilePath = OUT_PATH + "/stockholm";

        // Command to run RNALalifold
        List<String> command;

        // If input is a FASTA file, use --id-start to specify blockNumber
        if (inputFilePath.endsWith(".fasta")) {
            String regex = "block_(\\d+)_realigned";

            // Compile the pattern
            Pattern pattern = Pattern.compile(regex);

            // Create a matcher to find the pattern in the given string
            Matcher matcher = pattern.matcher(inputFilePath);
            // Check if the pattern matches and extract the block number
            int blockNumber = 0;
            if (matcher.find()) {
                // Group 1 contains the block number
                blockNumber = Integer.parseInt(matcher.group(1));
            } else {
                // Handle case where block number is not found
                throw new IllegalArgumentException("Block number not found in the input file path");
            }
            command = Arrays.asList(
                    ALIFOLDBINARY,
                    "--id-prefix=alifold",    // Prefix for output file names
                    "--id-start=" + blockNumber, // Start the ID at the specified blockNumber
                    "--noLP",                 // No lonely pairs
                    "--maxBPspan=300",        // Maximum base pair span
                    "--ribosum_scoring",      // Use Ribosum scoring
                    "--aln-stk",              // Output alignment in Stockholm format
                    inputFilePath             // Input FASTA file path
            );
        } else {
            // Command for MAF files (without --id-start since RNALalifold handles blocks automatically)
            command = Arrays.asList(
                    ALIFOLDBINARY,
                    "--id-prefix=alifold",    // Prefix for output file names
                    "--noLP",                 // No lonely pairs
                    "--maxBPspan=300",        // Maximum base pair span
                    "--ribosum_scoring",      // Use Ribosum scoring
                    "--aln-stk",              // Output alignment in Stockholm format
                    inputFilePath             // Input MAF file path
            );
        }

        ProcessBuilder pb = new ProcessBuilder(command);
        pb.directory(new File(outputFilePath));

        // Start the RNALalifold process
        Process process = pb.start();
        int exitCode = process.waitFor();

        if (exitCode != 0) {
            throw new IOException("RNALalifold exited with error code: " + exitCode);
        }

        if (VERBOSE && inputFilePath.endsWith(".fasta")) {
            System.out.println("RNALalifold output saved to: " + outputFilePath);
        }
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
}
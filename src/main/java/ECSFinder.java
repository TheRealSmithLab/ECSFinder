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
                    "  |  _| _ | |        _.____.    |  _|  [  | [ .-. |/ /'\\' |/ /__\\\\[ /'\\] \n" +
                    " _| |__/ |\\ .___.'\\| \\____) |  _| |_    | |  | | | || \\__/  || \\__., | |     \n" +
                    "|________| .____ .' \\______.' |_____|  [___][___||__]'.__.;__]'.__.'[___]    \n" +
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
        String nameOutputECS= "predicted_ECS.csv";
        mergeLogFiles(OUT_PATH + "/csv/", OUT_PATH + "/merged_hits.csv");
        callRScript(OUT_PATH + "/merged_hits.csv", OUT_PATH + "/"+nameOutputECS);

        File outputFile = new File(OUT_PATH + "/"+nameOutputECS);

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
                        File fileToDelete = new File(OUT_PATH+"/aln", nameFileWithUnderscores + ".aln");

                        // Check if the file exists and delete it
                        if (fileToDelete.exists()) {
                            fileToDelete.delete();
                        }
                    }

                    if (new File(OUT_PATH+"/aln").isDirectory()) {
                        if (new File(OUT_PATH+"/aln").list() != null && new File(OUT_PATH+"/aln").list().length == 0) {
                            System.out.println("No predicted ECS.");
                            new File(OUT_PATH+"/aln").delete();

                        }
                    }


                }

                // Close the reader
                reader.close();
            } catch (IOException e) {
                e.printStackTrace();
            }

            // Also delete the structure_input_sense.csv file
            File CsvFile = new File(OUT_PATH, "csv");

            // Check if the CSV file exists and delete it
            if (CsvFile.exists()) {
                CsvFile.delete();
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

    public static void processWithMafft(String inputFile) throws IOException, InterruptedException {
        // Convert MAF to FASTA and realign with MAFFT
        convertMafToSeparateFastas(inputFile, OUT_PATH);  // Converts to smaller blocks with overlap

        // Realign each block individually to avoid memory issues
        File outputDir = new File(OUT_PATH + "/output_fasta_dir");
        File[] fastaFiles = outputDir.listFiles((dir, name) -> name.endsWith(".fasta"));

        if (fastaFiles != null) {
            // Iterate through each FASTA block file
            for (File fastaFile : fastaFiles) {
                if(VERBOSE) {
                    System.out.println("Realigning file: " + fastaFile.getName());
                }
                // Realign each small block and save the result
                File realignedOutput = new File(fastaFile.getAbsolutePath().replace(".fasta", "_realigned.fasta"));
                realignSequences(fastaFile, outputDir);  // Pass block files to realignment
                runRNALalifold(realignedOutput.getAbsolutePath());  // Run RNA folding on the realigned blocks
            }
        }
      //  cleanUpFolder(outputDir);  // Optionally remove temporary files after processing
    }



    private static void processAlignmentBlocks(String inputFile) throws IOException, InterruptedException, ExecutionException {
        BufferedReader readFile = new BufferedReader(new FileReader(inputFile));
        int blockAln = 0;
        StringBuilder temp = new StringBuilder();
        ExecutorService multiThreads = Executors.newFixedThreadPool(NTHREDS);
        List<Future<?>> futures = new ArrayList<>();
        String path_aln = OUT_PATH + "/aln/";
        createDirectory(path_aln);
        String path_csv = OUT_PATH + "/csv/";
        createDirectory(path_csv);
        String line;
        while ((line = readFile.readLine()) != null) {
            if (line.length() >= 1 && line.charAt(0) != '#') {
                if (line.charAt(0) == 'a') {
                    blockAln++;
                } else if (line.charAt(0) == 's') {
                    temp.append(line).append("@");
                }
            }


            // Ensure temp has sufficient data before processing
            else if (!temp.toString().isEmpty() && (temp.toString().split("@").length >= 3) && line.equals("")) {
                processAlignmentBlock(temp.toString(), blockAln, futures, multiThreads);
                temp = new StringBuilder();  // Reset for the next block
            } else {
                temp = new StringBuilder();  // Reset if not enough sequences
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
        File dir = new File(OUT_PATH + "/stockholm");
        File[] stkFiles = dir.listFiles((dir1, name) -> name.endsWith(finalName + ".stk"));

        // Check if there are any Stockholm files found
        if (stkFiles == null || stkFiles.length == 0) {
            System.out.println("No Stockholm files found for alignment block " + blockAln + ". Not enough different species found in the alignments.");
            System.exit(1);  // Exit the program with a non-zero status
        }
        if(MAFFT) {
            // Regex pattern to match filenames with exactly two underscores
            String regex = "alifold_(\\d+)_" + finalName + "\\.stk";  // This will capture the number between the underscores

            // Compile the pattern
            Pattern pattern = Pattern.compile(regex);

            // Map to store file and the extracted number for sorting
            Map<File, Integer> fileNumberMap = new HashMap<>();

            // Iterate through all files
            for (File file : stkFiles) {
                Matcher matcher = pattern.matcher(file.getName());

                // Check if the file name matches the pattern
                if (matcher.matches()) {
                    // Extract the number between the underscores and convert it to an integer
                    int number = Integer.parseInt(matcher.group(1));
                    fileNumberMap.put(file, number);  // Store the file and the extracted number
                }
            }

            // Check if no matching files were found after filtering
            if (fileNumberMap.isEmpty()) {
                System.out.println("No matching Stockholm files found for alignment block " + blockAln + ". Not enough different species found in the alignments.");
                System.exit(1);  // Exit the program with a non-zero status
            }

            // Sort the files by the extracted number
            List<File> sortedFiles = new ArrayList<>(fileNumberMap.keySet());
            sortedFiles.sort(Comparator.comparingInt(fileNumberMap::get));  // Sort based on the extracted number

            // Iterate through the sorted files
            for (File file : sortedFiles) {
                if (VERBOSE) {
                    System.out.println("Processing file: " + file.getName());
                }
                processStockholmFile(file, mafTabTemp, associativeList, futures, multiThreads);
            }
        } else{
            for (File file : stkFiles) {
                processStockholmFile(file, mafTabTemp, associativeList, futures, multiThreads);
            }

            }
    }


    private static void processStockholmFile(File stockholmFile, String[] mafTabTemp, ArrayList<String[]> associativeList,
                                             List<Future<?>> futures, ExecutorService multiThreads) {
        String result = "";
        if(MAFFT) {
            // Regular expression to capture the desired parts
            String regex = "alifold_(\\d+)_(\\d+)";
            Pattern pattern = Pattern.compile(regex);

            // Match the string against the pattern
            Matcher matcher = pattern.matcher(stockholmFile.getName());

            int blockStart = 0;

            // Parse the filename and extract the first and second parts
            if (matcher.find()) {

                String firstPart = matcher.group(1);
                // Extract the second part and strip leading zeros
                String secondPart = matcher.group(2).replaceFirst("^0+", "");

                result = secondPart + "_" + firstPart;

                // Get the blockStart from the blockStartMap
                blockStart = blockStartMap.getOrDefault(result, 0); // Default to 0 if not found

            } else {
                System.out.println("Error: Filename does not match expected pattern: " + stockholmFile.getName());
                return; // Exit the method early if the filename does not match the expected format
            }
        } else {
            // Regular expression to capture the desired parts
            String regex = "alifold_(\\d+)";
            Pattern pattern = Pattern.compile(regex);

            // Match the string against the pattern
            Matcher matcher = pattern.matcher(stockholmFile.getName());

            int blockStart = 0;

            // Parse the filename and extract the first and second parts
            if (matcher.find()) {

                // Extract the second part and strip leading zeros
                String secondPart = matcher.group(1).replaceFirst("^0+", "");

                result = secondPart;


            }
        }
            // Open and process the Stockholm file
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

                    if (!associativeList.isEmpty() && currentLine.startsWith("//")) {
                        // Adjust motif coordinates using blockStart and pass it to processMotif
                        processMotif(mafTabTemp, arrayName, associativeList, gcReference, gcSScons, futures, multiThreads, result);
                    }
                }
            } catch (IOException e) {
                e.printStackTrace();
            }

    }



    private static void processMotif(String[] mafTabTemp, String[] arrayName, ArrayList<String[]> associativeList,
                                     String gcReference, String gcSScons, List<Future<?>> futures,
                                     ExecutorService multiThreads, String result) {
        // Adjust the motif start position by adding blockStart to ensure correct coordinates
        int[] cordMotif;
        if (MAFFT) {
           cordMotif = getRealCoordinates(Integer.parseInt(arrayName[4]), mafTabTemp, associativeList.get(0)[1], result);
        } else {
            cordMotif = getRealCoordinates(Integer.parseInt(arrayName[3]), mafTabTemp, associativeList.get(0)[1], result);
        }

        String loci = Arrays.toString(cordMotif);
        String chrom = mafTabTemp[1].substring((mafTabTemp[1].lastIndexOf(".") + 1));
        String lociChrm ="";
        if (MAFFT) {
            // Add blockStart to adjust loci positions properly
             lociChrm = chrom + ", " + loci.substring(1, loci.length() - 1) + ", " + mafTabTemp[4] + ", "
                    + arrayName[4] + ", " + arrayName[5] + ", " + gcReference + ", " + gcSScons;
        } else {
             lociChrm = chrom + ", " + loci.substring(1, loci.length() - 1) + ", " + mafTabTemp[4] + ", "
                    + arrayName[3] + ", " + arrayName[4] + ", " + gcReference + ", " + gcSScons;
        }

        String[] arrayLociChrm = lociChrm.split(", ");

        // Skip short alignments
        if (Integer.parseInt(arrayLociChrm[2]) - Integer.parseInt(arrayLociChrm[1]) < 50) {
            return;
        }

        // Create the ScanItFast task and submit it to the thread pool
        ScanItFast aln = new ScanItFast(associativeList, arrayLociChrm, new File(OUT_PATH), SSZBINARY, VERBOSE);
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
        if (line == null || line.trim().isEmpty()) {
            return null;  // Handle invalid lines
        }
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
     * Converts a MAF file to separate FASTA files, keeping all species aligned and splitting sequences into blocks of a defined length.
     *
     * @param mafFilePath   Path to the input MAF file.
     * @param outputDirPath Directory to output separate FASTA files.
     */
    public static void convertMafToSeparateFastas(String mafFilePath, String outputDirPath) {
        String fastaOutput = "output_fasta_dir";  // Directory where the FASTA files will be saved
        int blockCount = 0;  // To track the block number
        int overlapLength = 299; // Set overlap length to 299 bases
        int maxBlockSize = 5000;  // Define the maximum block size for each sequence

        try (BufferedReader reader = new BufferedReader(new FileReader(mafFilePath))) {
            // Create the directory if it does not exist
            Files.createDirectories(Paths.get(outputDirPath + "/" + fastaOutput));

            Map<String, StringBuilder> speciesSequences = new LinkedHashMap<>();  // Store species' sequences together
            String line;
            int currentBlockLength = 0;

            while ((line = reader.readLine()) != null) {
                line = line.trim();  // Trim whitespace

                if (line.startsWith("#") || line.isEmpty()) {
                    // Skip comments and empty lines
                    continue;
                }

                if (line.startsWith("a")) {
                    // When starting a new alignment block, reset and write the previous block if it exists
                    if (!speciesSequences.isEmpty()) {
                        blockCount++;
                        splitAndWriteBlocks(outputDirPath, blockCount, speciesSequences, maxBlockSize, overlapLength);
                        speciesSequences.clear();  // Clear for the next block
                    }
                    currentBlockLength = 0;  // Reset block length
                } else if (line.startsWith("s")) {
                    // Process sequence lines
                    String[] tokens = line.split("\\s+");
                    String sequenceId = tokens[1];  // Sequence ID
                    String sequence = tokens[tokens.length - 1];  // Aligned sequence with gaps

                    // Add or append the sequence to the current block for the species
                    speciesSequences.computeIfAbsent(sequenceId, k -> new StringBuilder()).append(sequence);
                    currentBlockLength = Math.max(currentBlockLength, speciesSequences.get(sequenceId).length());
                }
            }

            // Write the last block if any sequences remain
            if (!speciesSequences.isEmpty()) {
                blockCount++;
                splitAndWriteBlocks(outputDirPath, blockCount, speciesSequences, maxBlockSize, overlapLength);
            }

        } catch (IOException e) {
            System.err.println("Error processing the MAF file: " + e.getMessage());
        }
    }

    /**
     * Splits the sequences in each species group into blocks of defined length and writes them to separate FASTA files.
     *
     * @param outputDirPath     Directory where the output FASTA files will be saved.
     * @param blockCount        The current block number (used for naming the files).
     * @param speciesSequences  Map of species to their sequences in the alignment block.
     * @param maxBlockSize      The maximum allowed length of each block.
     * @param overlapLength     The number of bases to overlap between blocks.
     * @throws IOException If an I/O error occurs.
     */
    // Map to store block start and part number information for each generated file
    private static Map<String, Integer> blockStartMap = new HashMap<>();


    private static void splitAndWriteBlocks(String outputDirPath, int blockCount, Map<String, StringBuilder> speciesSequences,
                                            int maxBlockSize, int overlapLength) throws IOException {

        int blockPart = 1;  // Keep track of block part (if we split it)
        boolean isSplit = false;

        // Get the length of sequences (they should all have the same length)
        int fullLength = speciesSequences.values().iterator().next().length();
        int blockStart = 0;  // This will track the starting position of each part within the block

        for (int i = 0; i < fullLength; i += (maxBlockSize - overlapLength)) {
            String blockFileName = outputDirPath + "/output_fasta_dir/block_" + blockCount + "_part_" + blockPart + ".fasta";
            blockStart = i;  // Track the start position for this part

            File blockFile = new File(blockFileName);

            // Store the block start and block part in the maps for later use
            blockStartMap.put(blockCount+"_"+blockPart, blockStart);

            try (BufferedWriter writer = new BufferedWriter(new FileWriter(blockFileName))) {
                for (Map.Entry<String, StringBuilder> entry : speciesSequences.entrySet()) {
                    String speciesId = entry.getKey();
                    String sequence = entry.getValue().toString();

                    // Extract the sub-sequence for this part
                    int end = Math.min(i + maxBlockSize, sequence.length());
                    String subSequence = sequence.substring(i, end);

                    // Write the sub-sequence to the FASTA file
                    writer.write(">" + speciesId + "\n");
                    writer.write(subSequence + "\n");
                }
            }
            if(VERBOSE) {
                System.out.println("Wrote block part: " + blockFileName + " starting at position: " + blockStart);
            }
            blockPart++;  // Increment the block part for the next segment

            if (i + maxBlockSize < fullLength) {
                isSplit = true;
            }
        }

        if (isSplit && VERBOSE) {
            System.out.println("Block " + blockCount + " was split into " + (blockPart - 1) + " parts.");
        }
    }


    /**
     * Realigns sequences using MAFFT and writes the output directly to a FASTA file.
     * @param inputFilePath    The path to the input FASTA file.
     * @param outputFilePath   The path to the output FASTA file.
     * @throws IOException           If an I/O error occurs.
     * @throws InterruptedException  If the MAFFT process is interrupted.
     */
    private static Map<String, String> homoSapiensSequences = new HashMap<>();
    public static void realignSequences(File inputFilePath, File outputFilePath) throws IOException, InterruptedException {
        String realignedFilePath = inputFilePath.getAbsolutePath().replace(".fasta", "_realigned.fasta");
        // Regular expression to extract the block number and part number
        String regex = "block_(\\d+)_part_(\\d+)_realigned\\.fasta";
        Pattern pattern = Pattern.compile(regex);
        Matcher matcher = pattern.matcher(realignedFilePath);
        String result="";
        if (matcher.find()) {
            // Extract the block number and part number
            String blockNumber = matcher.group(1);
            String partNumber = matcher.group(2);

            // Combine them in the desired format
             result = blockNumber + "_" + partNumber;

        } else {
            System.out.println("No match found.");
        }

        // Run MAFFT on the gap-stripped sequences
        List<String> command = Arrays.asList(MAFFTBINARY, "--quiet", "--thread", String.valueOf(NTHREDS), inputFilePath.getAbsolutePath());
        ProcessBuilder pb = new ProcessBuilder(command);
        Process mafftProcess = pb.start();

        // Capture the realigned sequences
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

                        // If this is the Homo sapiens sequence, update the map with the realigned sequence
                        if (speciesReal.toLowerCase().contains("homo")) {
                            homoSapiensSequences.put(result, sequence.toString());
                        }

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
        File tempInputFile = new File(inputFilePath.getParent(), "temp_gap_stripped.fasta");
        if (tempInputFile.exists()) {
            tempInputFile.delete();
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
            String regex = "block_(\\d+)_part_(\\d+)_realigned";

            // Compile the pattern
            Pattern pattern = Pattern.compile(regex);

            // Create a matcher to find the pattern in the given string
            Matcher matcher = pattern.matcher(inputFilePath);

            // Variables to store block number and part number
            int blockNumber = 0;
            int partNumber = 0;

            // Check if the pattern matches and extract the block and part numbers
            if (matcher.find()) {
                blockNumber = Integer.parseInt(matcher.group(1));  // Extract and convert the block number
                partNumber = Integer.parseInt(matcher.group(2));
            } else {
                // Handle case where block number is not found
                throw new IllegalArgumentException("Block number not found in the input file path");
            }
            command = Arrays.asList(
                    ALIFOLDBINARY,
                    "--id-prefix=alifold_"+String.valueOf(partNumber),    // Prefix for output file names
                    "--id-start=" + String.valueOf(blockNumber), // Start the ID at the specified blockNumber
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

        // Handle stdout and stderr
        BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
        BufferedReader errorReader = new BufferedReader(new InputStreamReader(process.getErrorStream()));

        // Create threads to handle both streams
        Thread outputThread = new Thread(() -> {
            // Read stdout but do not print it
            reader.lines().forEach(line -> {

            });
        });

        Thread errorThread = new Thread(() -> {
            errorReader.lines().forEach(line -> {

            });
        });


        // Start both threads
        outputThread.start();
        errorThread.start();
        int exitCode = process.waitFor();

        // Ensure the threads have finished reading
        outputThread.join();
        errorThread.join();
        if (exitCode != 0) {
            throw new IOException("RNALalifold exited with error code: " + exitCode);
        }

        if (VERBOSE && inputFilePath.endsWith(".fasta")) {
            System.out.println("RNALalifold output saved to: " + outputFilePath);
        }
    }

    private static int[] getRealCoordinates(int start, String[] mafCord, String motifHuman, String blockPartKey) {
        int[] cordFinal = new int[2];
        int[] cordFinalPlus1 = new int[2];
        int nuc=0;
        int nuclStockholm=0;
        if(MAFFT) {
            // Retrieve the Homo sapiens realigned sequence and the block start
            String homoSapiensSequence = homoSapiensSequences.get(blockPartKey);
            String homoRealigned = homoSapiensSequence.substring(0, start);
             nuclStockholm = homoRealigned.replaceAll("-", "").length();  // Nucleotide count up to 'start'

            int blockStart = blockStartMap.getOrDefault(blockPartKey, 0);
            // Get the corresponding original MAF sequence (ignoring gaps)
            String withoutGap = mafCord[6].substring(0, blockStart);
             nuc = withoutGap.replaceAll("-", "").length();  // Count nucleotides without gaps from the original sequence
        } else{
            String withoutGap = mafCord[6].substring(0, start);
            nuc = withoutGap.replaceAll("-", "").length();
            nuclStockholm=0;
        }
        // Forward strand case
        if (mafCord[4].equals("+")) {
            int lociStart = Integer.parseInt(mafCord[2]) + nuc+nuclStockholm;  // Adjust with original MAF nucleotide count
            int lociEnd = lociStart + motifHuman.replaceAll("-", "").length();  // Add motif length excluding gaps
            cordFinal = new int[]{lociStart, lociEnd};
        }
        // Reverse strand case
        else {
            int lociEnd = (Integer.parseInt(mafCord[5])+1 - (Integer.parseInt(mafCord[2]) + nuc+nuclStockholm)) + 1;  // Adjust end coordinate
            int lociStart = lociEnd - motifHuman.replaceAll("-", "").length();  // Subtract motif length
            cordFinal = new int[]{lociStart, lociEnd};
        }

        cordFinalPlus1[0] = cordFinal[0];
        cordFinalPlus1[1] = cordFinal[1] - 1;  // Adjust the end for 0-based coordinates
        return cordFinalPlus1;
    }

    public static void mergeLogFiles(String logDirPath, String finalCsvPath) {
        File logDir = new File(logDirPath);
        File[] logFiles = logDir.listFiles((dir, name) ->  name.endsWith(".csv"));

        if (logFiles != null && logFiles.length > 0) {
            try (BufferedWriter finalWriter = new BufferedWriter(new FileWriter(finalCsvPath))) {
                finalWriter.write("name_file,min_energy,pseudo_energy,log_min_evalue,covarying_bp,MPI,average_MFE_sample,sd_sample,zscore\n");

                for (File logFile : logFiles) {
                    try (BufferedReader logReader = new BufferedReader(new FileReader(logFile))) {
                        String line;
                        while ((line = logReader.readLine()) != null) {
                            finalWriter.write(line);
                            finalWriter.newLine();
                        }
                    } catch (IOException e) {
                        System.err.println("Error reading log file: " + logFile.getName());
                        e.printStackTrace();
                    }
                    // Delete the log file after merging
                    logFile.delete();
                }
            } catch (IOException e) {
                System.err.println("Error writing to final CSV file.");
                e.printStackTrace();
            }
        } else {
            System.out.println("No log files found to merge.");
        }
    }


}
import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

public class DataLoader {

    public final File[] rfamFolder;
    public Map<String, List<String[]>> sequenceNames;
    public List<String> totalSequence;
    public Map<String, List<String[]>> sequenceChars;
    public List<boolean[]> beenSampled;
    private List<boolean[]> isWrongSize;
    public List<String> totalRfam;
    private File mafFolder;
    private File mafSourceFolder;
    private int minPerFam;


    public DataLoader(File mafFolder, File[] rfamFolder, File mafSourceFolder, int minPerFam) {
        this.mafFolder = mafFolder;
        this.rfamFolder = rfamFolder;
        this.mafSourceFolder = mafSourceFolder;
        this.minPerFam = minPerFam;
        this.sequenceChars = new HashMap<>();
        this.beenSampled = new ArrayList<>();
        this.isWrongSize = new ArrayList<>();
        this.totalSequence = new ArrayList<>();
        this.totalRfam = new ArrayList<>();
        this.sequenceNames = new HashMap<>();
    }
/**
 * Loads data from input files and processes them according to specified criteria.
 * @param files Array of Files to be processed.
 * @param verbose Flag to print detailed process information.
 */


    /**
     * Processes FASTA files, reading sequences and storing relevant information.
     * @param files Array of input FASTA files.
     * @param verbose Flag for verbose output.
     */
    public void processFastaFiles(File[] files, boolean verbose) throws IOException {
        for (File file : files) {
            if (verbose) {
                System.err.println("Parsing file: " + file.getName());
            }


            Map<String, String> rfamAln = new TreeMap<>();
            try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
                String line;
                String sequence = "";
                String species = "";
                boolean isFirstSequence = true;

                while ((line = reader.readLine()) != null) {
                    if (line.startsWith(">")) {

                        if (!isFirstSequence) {

                            if (sequence.replaceAll("-", "").length() > 50) {
                                rfamAln.putIfAbsent(species, sequence.toUpperCase());
                            }
                        } else {
                            isFirstSequence = false;
                        }

                        species = line.substring(1).trim();
                        sequence = "";
                    } else {
                        sequence += line.trim();
                    }
                }

                if (!isFirstSequence && sequence.replaceAll("-", "").length() > 50) {
                    rfamAln.put(species, sequence.toUpperCase());
                }
            } catch (IOException e) {
                System.err.println("Failed to read file: " + file.getName());
                throw e;
            }


            storeSequenceData(rfamAln, verbose, file.getName());
        }
    }


    /**
     * Stores sequence data into class members from a given RFAM alignment map.
     * @param rfamAln Map containing species as keys and sequences as values.
     * @param verbose Flag for verbose output.
     */
    private void storeSequenceData(Map<String, String> rfamAln, boolean verbose, String fileName) {
        if (!rfamAln.isEmpty()) {

            String rfamId = fileName.substring(0, fileName.lastIndexOf(".fasta"));// You need to implement this method
            List<String[]> sequenceNamePairs = rfamAln.entrySet().stream()
                    .map(entry -> new String[]{entry.getValue(), entry.getKey()})
                    .collect(Collectors.toList());

            // Store the sequence-name pairs in sequenceNames under the rfamId
            sequenceNames.put(rfamId, sequenceNamePairs);


            totalRfam.add(rfamId);
            totalSequence.addAll(rfamAln.values());

            // Initialize beenSampled and isWrongSize if needed
            boolean[] sampledFlags = new boolean[sequenceNamePairs.size()];
            Arrays.fill(sampledFlags, false);
            beenSampled.add(sampledFlags);

            boolean[] sizeFlags = new boolean[sequenceNamePairs.size()];
            Arrays.fill(sizeFlags, false);
            isWrongSize.add(sizeFlags);

            if (verbose) {
                System.err.println("Loaded " + sequenceNamePairs.size() + " sequences for RFAM ID: " + rfamId);
            }
        }
    }

    public void processMafFiles(Map<String, List<String[]>> sampledSequences, File outFile, boolean verbose) {
        File targetDir = new File(mafFolder.getPath());
        if (!targetDir.exists() || !targetDir.isDirectory()) {
            System.err.println("Target directory does not exist or is not a directory: " + targetDir.getPath());
            return;
        }


        File[] filesToProcess = mafSourceFolder.listFiles((dir, name) -> {
            if (name.endsWith(".fasta")) {
                int lastUnderscore = name.lastIndexOf('_');
                int secondLastUnderscore = name.lastIndexOf('_', lastUnderscore - 1);
                if (secondLastUnderscore != -1 && lastUnderscore != -1) {
                    try {
                        String numberStr = name.substring(secondLastUnderscore + 1, lastUnderscore);
                        int number = Integer.parseInt(numberStr);
                        return number >= minPerFam;
                    } catch (NumberFormatException e) {
                        return false;
                    }
                }
            }
            return false;
        });

        if (filesToProcess == null || filesToProcess.length == 0) {
            if (verbose) {
                System.out.println("No files found in source directory for processing.");
            }
            return;
        }

        List<String> keysList = new ArrayList<>(sampledSequences.keySet());

        for (int i = 0; i < filesToProcess.length && i < keysList.size(); i++) {
            File currentFile = filesToProcess[i];
            String key = keysList.get(i);
            List<String[]> sequencesToProcess = sampledSequences.getOrDefault(key, Collections.emptyList());

            if (!sequencesToProcess.isEmpty() && verbose) {
                System.out.println("Processing " + sequencesToProcess.size() + " sequences for file " + currentFile.getName() + " using key: " + key);
            }

            try {
                processFile(sequencesToProcess, currentFile, verbose, key, outFile);

                File newLocation = new File(targetDir, currentFile.getName());
                if (currentFile.renameTo(newLocation)) {
                    if (verbose) {
                        System.out.println("Successfully moved file: " + currentFile.getName());
                    }
                } else {
                    if (verbose) {
                        System.err.println("Failed to move file: " + currentFile.getName());
                    }
                }
            } catch (IOException e) {
                if (verbose) {
                    System.err.println("Failed to process file: " + currentFile.getName() + " due to an error: " + e.getMessage());
                }
            }
        }
    }



    private void processFile(List<String[]> sampledSequencesWithNames, File file, boolean verbose, String key, File outFile) throws IOException {
        int number = 0;
        List<String> originalLines = Files.readAllLines(file.toPath());
        List<String[]> modifiedSequencesWithNamesGapless = new ArrayList<>();
        List<String[]> modifiedSequencesWithNamesModified = new ArrayList<>();// To store both sequences and names
        List<String> mafSpecies =new ArrayList<>();
        StringBuilder sequenceBuilder = new StringBuilder();
        String currentHeader = null;
        String currentName = ""; // To store the current sequence's name
        File positionsFile = new File(outFile+"/"+"positions.txt");
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(positionsFile, true))) { // Append mode
            for (String line : originalLines) {
                if (line.startsWith(">")) {
                    // Extract the name from the header line, assuming it follows ">" directly
                    currentName = line.substring(1).trim();

                    // Process the previous sequence
                    if (currentHeader != null && sequenceBuilder.length() > 0 && number - 1 < sampledSequencesWithNames.size()) {
                        mafSpecies.add(currentHeader);
                        String[] insertResult = insertAndCutSequence(sequenceBuilder.toString(), sampledSequencesWithNames.get(number - 1)[0], 300);
                        String modifiedSequence = insertResult[0];
                        String gaplessSequence = modifiedSequence.replace("-", "").replace("U", "T");
                        // Writing positions and key to the file
                        writer.write("Key: " + key + ", Start position: " + insertResult[1] + ", Cut position: " + insertResult[2]);
                        writer.newLine();
                        // Add both the modified sequence and its name
                        modifiedSequencesWithNamesGapless.add(new String[]{gaplessSequence, sampledSequencesWithNames.get(number - 1)[1], currentHeader});
                        modifiedSequencesWithNamesModified.add(new String[]{modifiedSequence, sampledSequencesWithNames.get(number - 1)[1], currentHeader});
                    }
                    number++;
                    currentHeader = line.trim();
                    sequenceBuilder = new StringBuilder();
                } else {
                    sequenceBuilder.append(line.trim());
                }
            }
            // Process the last sequence if there's one
            if (currentHeader != null && sequenceBuilder.length() > 0 && number - 1 < sampledSequencesWithNames.size()) {
                String[] insertResult = insertAndCutSequence(sequenceBuilder.toString(), sampledSequencesWithNames.get(number - 1)[0], 300);
                String modifiedSequence = insertResult[0];
                modifiedSequencesWithNamesModified.add(new String[]{modifiedSequence, sampledSequencesWithNames.get(number - 1)[1], currentHeader});
            }

            writeSequencesToFile(modifiedSequencesWithNamesModified, file, key, verbose);
        }

        // return modifiedSequencesWithNames;
    }

    private void writeSequencesToFile( List<String[]> modifiedSequencesWithNamesModified,File originalFile, String key, boolean verbose) throws IOException {
        // Filenames based on the original file name and the key

        Path modifiedPath =Paths.get(mafFolder+ "/preMafft");;
        Path gaplessPath = Paths.get(mafFolder+"/gapless");

        // Assuming originalFile is a File object
        String originalFileName = originalFile.getName();
        String modifiedFileName = "modified_" +key+"_"+ originalFileName; // Example of modifying the file name
        Path modifiedFilePath = modifiedPath.resolve(modifiedFileName);
        Path gaplessFilePath = gaplessPath.resolve("gapless_" +key+"_"+ originalFileName);

        Files.createDirectories(modifiedPath);
        Files.createDirectories(gaplessPath);

        try (BufferedWriter modifiedWriter = Files.newBufferedWriter(modifiedFilePath);
             BufferedWriter gaplessWriter = Files.newBufferedWriter(gaplessFilePath)) {

            for (String[] sequenceWithName : modifiedSequencesWithNamesModified) {
                String species= sequenceWithName[2].split("\\.")[0];
                String numb=sequenceWithName[2].split("\\.")[1];
                String name = species+sequenceWithName[1].split("[./_]", 2)[0]+"."+numb;
                String sequence = sequenceWithName[0]; // Sequence at index 0
                String gaplessSequence = sequence.replace("-", "").replace("U", "T"); // Process for gapless

                // Write to modified sequences file
                modifiedWriter.write( name+"_"+key + "\n");
                modifiedWriter.write(sequence + "\n\n");

                // Write to gapless sequences file
                gaplessWriter.write(name+"_"+key + "\n");
                gaplessWriter.write(gaplessSequence + "\n\n");
            }
        }

        if (verbose) {
            System.out.println("Sequences written to:\n- " + modifiedPath + "\n- " + gaplessPath);
        }
    }





    private String[] insertAndCutSequence(String originalSequence, String insertSequence, int position) {
        int realPosition = calculatePositionIgnoringGaps(originalSequence, position);
        String start = originalSequence.substring(0, realPosition);
        String end = originalSequence.substring(realPosition);
        String modified = start + insertSequence + end;
        int cutPosition = calculatePositionIgnoringGaps(modified, position + insertSequence.length() + 300);
        // Returning the modified sequence along with start and end positions
        return new String[]{modified.substring(0, cutPosition), String.valueOf(realPosition), String.valueOf(cutPosition)};
    }


    private int calculatePositionIgnoringGaps(String sequence, int n) {
        int count = 0;
        for (int i = 0; i < sequence.length(); i++) {
            if (sequence.charAt(i) != '-') {
                count++;
                if (count == n) return i + 1;
            }
        }
        return sequence.length();
    }

}

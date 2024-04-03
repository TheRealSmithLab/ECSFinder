import java.io.IOException;
import java.util.*;

public class SequenceSampler {
    private static boolean verbose;
    private int minPerFam;
    private int maxPerFam;
    private int totalSeqs;
    private int MIN_PI;
    private int MAX_PI;
    private Random random = new Random();
    public Map<String, List<String[]>> sequenceNames;

    public SequenceSampler(int minPerFam, int maxPerFam, int minPi, int maxPi, int totalSeqs, boolean verbose) {
        this.minPerFam = minPerFam;
        this.maxPerFam = maxPerFam;
        this.MIN_PI = minPi;
        this.MAX_PI = maxPi;
        this.totalSeqs = totalSeqs;
        this.verbose = verbose;
        this.sequenceNames=sequenceNames;
    }

    public Map<String, List<String[]>> sampleSequences(Map<String, List<String[]>> sequenceNames) throws IOException {
        Map<String, List<String[]>> sampledFamilies = new HashMap<>();
        int sequencesSampled = 0;

        List<String> shuffledFamilies = new ArrayList<>(sequenceNames.keySet());
        Collections.shuffle(shuffledFamilies, random); // Shuffle families for randomized access

        for (String family : shuffledFamilies) {
            if (sequencesSampled >= totalSeqs) break; // Stop if total sequence count is reached

            List<String[]> familySequences = new ArrayList<>(sequenceNames.get(family));
            Collections.shuffle(familySequences, random); // Shuffle sequences within a family

            if (familySequences.size() < minPerFam) {
                if (verbose) System.out.println("Skipping family " + family + " due to insufficient sequences.");
                continue;
            }

            List<String[]> sampledForFamily = new ArrayList<>();

            // Iterate through the shuffled list of sequences for this family
            for (String[] sequenceWithName : familySequences) {
                if (sampledForFamily.size() >= maxPerFam || sequencesSampled >= totalSeqs) break; // Enforce limits

                String sequence = sequenceWithName[0];
                double pid = calculatePID(familySequences.get(0)[0], sequence); // Compute pairwise identity using the first sequence as reference
                sequencesSampled++;
                if ((pid >= MIN_PI ) && (pid <= MAX_PI+20)) {
                    sampledForFamily.add(sequenceWithName); // Add this sequence and name if within PID bounds

                }
            }

            // Only add sampledForFamily to sampledFamilies if it contains at least maxPerFam sequences
            if (sampledForFamily.size() >= maxPerFam) {
                sampledFamilies.put(family, sampledForFamily);

                if (verbose) {
                    System.out.println("Sampled " + sampledForFamily.size() + " sequences from family " + family + ".");
                }
            }
        }

        return sampledFamilies;
    }



    private double calculatePID(String seq1, String seq2) {
        int matches = 0, nonGapCount = 0;

        for (int i = 0; i < seq1.length(); i++) {
            if (seq1.charAt(i) != '-' && seq2.charAt(i) != '-') {
                nonGapCount++;
                if (seq1.charAt(i) == seq2.charAt(i)) {
                    matches++;
                }
            }
        }

        return nonGapCount > 0 ? (double) matches / nonGapCount * 100.0 : 0.0;
    }
}

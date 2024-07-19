import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.regex.*;
import java.util.stream.*;

public class FilterOutput {
    double min_eval;
    public double processFilesWithSuffix(String directoryPath, String suffix, String searchPattern) {

        try {
            // Get a list of all files ending with the given suffix in the specified directory
            List<Path> files = Files.list(Paths.get(directoryPath))
                    .filter(path -> path.toString().endsWith(suffix))
                    .collect(Collectors.toList());

            // Loop through each file
            for (Path filePath : files) {
                List<Double> eValues= new ArrayList<>();
                BufferedReader reader = new BufferedReader(new FileReader(filePath.toFile()));
                String line;

                // Extract values based on the search pattern
                while ((line = reader.readLine()) != null) {
                    Matcher matcher = Pattern.compile(searchPattern + "([^ ]*)").matcher(line);
                    if (matcher.find()) {
                        eValues.add(Double.valueOf(matcher.group(1)));
                    }
                }
                reader.close();
                min_eval=Collections.min(eValues);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return min_eval;
    }

    public double[] processTxtFiles(String directoryPath) {
        double[] energies = new double[2];
        try {
            // Get a list of all files ending with '.txt' in the specified directory
            List<Path> txtFiles = Files.list(Paths.get(directoryPath))
                    .filter(path -> path.toString().endsWith(".txt"))
                    .collect(Collectors.toList());

            // Loop through each file
            for (Path filePath : txtFiles) {
                BufferedReader br = new BufferedReader(new FileReader(filePath.toFile()));
                String line;
                Pattern pattern = Pattern.compile("\\(([^=]+)=([^+]+)\\+([^\\)]+)\\)");

                while ((line = br.readLine()) != null) {
                    Matcher matcher = pattern.matcher(line);
                    if (matcher.find()) {
                        String part2 = matcher.group(2).trim().replaceAll("\\)", "");
                        String part3 = matcher.group(3).trim().replaceAll("\\)", "");
                         energies = new double[]{Double.parseDouble(part2), Double.parseDouble(part3)};
                    }
                }
                br.close();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return energies;
    }
}

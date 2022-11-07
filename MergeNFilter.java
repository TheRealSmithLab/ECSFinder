import java.util.*;import java.io.*;
/**
 * This program filters multiple sequence alignments and remove duplicates species within alignment blocks
 * @author Vanda Gaonach-Lovejoy 19/03/21
 */

public class MergeNFilter {

	public static void main(String[] args) throws IOException {


		if (args.length == 0) {
			System.out.println("Warning: maf files are expected");
			System.exit(-1);
		}

		if (args[0].split("/").length==1){
			System.out.println("Please enter the full path of the input file");
			System.exit(0);
		}

		try {

			int input = 0;
			BufferedWriter segDups = new BufferedWriter(new FileWriter(args[0].substring(0, args[0].lastIndexOf("."))
					+ "-removedLines.txt"));
			BufferedWriter out = new BufferedWriter(new FileWriter(args[0].substring(0, args[0].lastIndexOf(
					".")) + "-output.maf"));
			;
			out.write("##maf version=2 \n" +
					"# original dump date: 2020-12-10 13:44:00\n# ensembl release: 103\n" +
					"# emf comment: Alignments: 46 eutherian mammals EPO\n# emf comment: Region:" +
					" Homo sapiens chromosome:GRCh38\n");

			while (input != args.length) {
				String file = args[input];
				BufferedReader entry = new BufferedReader(new FileReader(file));
				String line;
				while ((line = entry.readLine()) != null) {

					if(line.length() != 0 && line.charAt(0) == 's') {
						LinkedHashMap<String, String[]> speciesSequences = new LinkedHashMap<>();
						ArrayList<String[]> duplicateSequences = new ArrayList<>();
						while (line.length() != 0 && line.charAt(0) == 's') {
							String[] arraySequenceInfo = line.split("\\s+");
							String nameSpeciesWithChro = arraySequenceInfo[1];
							String nameSpeciesOnly = nameSpeciesWithChro.substring(0,
									nameSpeciesWithChro.indexOf("."));
							if (!nameSpeciesOnly.equals("ancestral_sequences")) {
								if (!(speciesSequences.containsKey(nameSpeciesOnly))) {
									speciesSequences.put(nameSpeciesOnly, arraySequenceInfo);
								} else {
									duplicateSequences.add(arraySequenceInfo);
								}
							}
							line = entry.readLine();
							if (line == null)
								System.exit(0);
						}
						if (speciesSequences.size() == 1) {
							duplicateSequences.add((speciesSequences.get("homo_sapiens")));
						}

						for (int i = 0; i < duplicateSequences.size(); i++) {
							String[] arraySpecies = duplicateSequences.get(i);
							if (i == 0) {
								segDups.write("\na score=0\n");
							}
							String stringSpecies = Arrays.toString(arraySpecies);
							String noBrackets = stringSpecies.replace("[", "")
									.replace("]", "")
									.replace(",", "\t");


							segDups.write(noBrackets + "\n");
						}

						if (speciesSequences.size() > 1){
							out.write( "a" +"\n");
							for (String key : speciesSequences.keySet()) {
								String[] value = speciesSequences.get(key);
								String eachSpeciesInfo = Arrays.toString(value);

								//remove the right and left bracket
								String noBrackets = eachSpeciesInfo.replace("[", "")
										.replace("]", "")
										.replace(",", "\t");
								out.write(noBrackets + "\n");
							}
						}
					}
					out.write("\n");
				}
				input++;
				entry.close();
			}
			out.write("\n");
			out.write("\n");
			out.close();
			segDups.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}

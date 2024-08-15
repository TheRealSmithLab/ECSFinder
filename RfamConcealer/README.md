# RfamConcealer

## Overview

`RfamConcealer` is a Java program that selects and processes random (sub)sequences from a reference database, specifically focusing on RFAM families. The tool can generate MAF files, realign sequences using MAFFT, and produce Clustal files for multiple sequence alignments.

## Features

- **Data Loading and Preparation**: Load and prepare sequences from FASTA files in the specified input directory.
- **Sequence Sampling**: Sample sequences from each RFAM family based on specified parameters.
- **MAFFT Realignment**: Realign sequences using MAFFT and save the results in MAF format.
- **Clustal File Generation**: Generate Clustal files for multiple sequence alignment.

## Requirements

- **Java 8** or higher
- **MAFFT** installed (with the `mafft-ginsi` executable available)
- A directory containing FASTA files as the input database

## Usage
You can run the RfamConcealer program from the command line using the following syntax: 'java RfamConcealer [options] -i [path to .fasta database folder]'
- **`-i <string>`**: (Required) Specifies the path to the directory containing the FASTA files.
- **`-v`**: Enables verbose output for detailed logging during execution.
- **`-min_pi <int>`**: Sets the minimum pairwise identity percentage (default: `10`).
- **`-max_pi <int>`**: Sets the maximum pairwise identity percentage (default: `80`).
- **`-max_f <int>`**: Sets the maximum number of sequences per RFAM family/input file (default: `20`).
- **`-min_f <int>`**: Sets the minimum number of sequences per RFAM family/input file (default: `2`).
- **`-t <int>`**: Sets the maximum total number of sequences to include (default: `200`).
- **`-o <string>`**: Specifies the output directory where the results will be saved (default: `./rfam_subset.fa`).

### Output

Upon successful execution, the program generates the following output:

- **MAF File**:
  - **Location**: The MAF file is saved in the output directory specified by the `-o` option (default: `./rfam_subset.fa/output.maf`).
  - **Description**: Contains realigned sequences in MAF format, including information on pairwise identities.

- **Clustal Files**:
  - **Location**: Stored in a `ClustalFiles` subdirectory within the specified output directory.
  - **Description**: Clustal format files generated for each RFAM family, useful for multiple sequence alignments.

- **Log Information**:
  - **Description**: If the `-v` (verbose) option is enabled, detailed logs about the processing steps are printed to the console, including sequence sampling and alignment information.

###Example

'-o RFAM_65MPI_80MPI_20species_new_2 -min_pi 65 -max_pi 80 -min_f 20 -max_f 20 -t 100000 -i /home/vandalovejoy/rfam_families_new/filtered_families_rfam/over20'
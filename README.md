 # ECSFinder

Scans multiple alignments for conserved RNA structures. Reads a set of maf files, calculates stats, scans with SISSIz or RNAz, outputs bed coordonates of high-confidence predictions


## Usage

```
java ECSFinder [options] -o output/directory -i input.maf (last parameter must be -i)
 Options:
   -c int number of CPUs for calculations (default 4)
   -g int max gap percentage of sequences for 2D prediction (default 50)
   -sszr double report SISSIz+RIBOSUM hits below this Z-score (default -3)
   -v verbose (messy but detailed) output
   -r R-scape option

```

 ## Output: Two types of results are produced:<br />
                            &ensp;&ensp;   (1)  the multiple sequence alignments associated to significant predictions<br />
                            &ensp;&ensp;        are saved to files in the folder specified with the -o option.<br />
                            &ensp;&ensp;        File names correspond to their genomic coordinates in a .bed-compatible format. <br />
                            &ensp;&ensp;        Ex: output/directory/chrX_12345_12500_80:75:23:14:8:z_300_+.aln<br />
                            &ensp;&ensp;   (2)  The genomic coordinates (.bed format) of ECSs are also written to the SDOUT<br />
                            &ensp;&ensp;        (see additional options below).<br />
 
 
 
 

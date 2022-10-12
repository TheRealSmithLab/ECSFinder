


<p align="center" width="100%">
    <img width="25%" src="https://user-images.githubusercontent.com/44384386/195381940-680064be-d53a-45b6-a5e1-a80ff1cb804e.jpg"> 
</p>



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

```

 ## Output: Two types of results are produced:
                              (1)  the multiple sequence alignments associated to significant predictions
                                   are saved to files in the folder specified with the -o option
                                    File names correspond to their genomic coordinates in a .bed-compatible format
                                    Ex: output/directory/chrX_12345_12500_80:75:23:14:8:z_300_+.aln
                              (2)  The genomic coordinates (.bed format) of ECSs are also written to the SDOUT
                                    (see additional options below)
 
 
 

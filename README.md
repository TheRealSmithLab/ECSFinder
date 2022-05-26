 ECSFinder

SCAN MULTIPLE ALIGNMENTS FOR CONSERVED RNA STRUCTURES
Reads a set of maf files, calculates stats, scans with SISSIz or RNAz, outputs bed coordonates of high-confidence predictions
                    Usage:     java  ECSFinder.jar [options] -o output/directory -i input.maf (last parameter must be -i)
                    Output: 	Two types of results are produced:
                               (1)  the multiple sequence alignments associated to significant predictions
                                    are saved to files in the folder specified with the -o option.
                                    File names correspond to their genomic coordinates in a .bed-compatible format. 
                                    Ex: output/directory/chrX_12345_12500_80:75:23:14:8:z_300_+.aln
                               (2)  The genomic coordinates (.bed format) of ECSs are also written to the SDOUT
                                    (see additional options below).
                    output: Number chromosome start loci  end loci number of species MPI standard deviation Normalized shannon entropy
                            GC % content Gap %content SISSIZ Z-score
                               ***  N.B. the score field corresponds to the SISSIz Z-score x-100
                    Options:
                      -all             write out bed entires for all sampled alignments to STDOUT
                      -c     int       number of CPUs for calculations (default 4)
                      -g     int       max gap percentage of sequences for 2D prediction (default 50)
                      -sszr  double    report SISSIz+RIBOSUM hits below this Z-score (default -3)
                      -v               verbose (messy but detailed) output
                      -r                R-scape to be used

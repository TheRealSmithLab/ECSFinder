 **ECSFinder

**SCAN MULTIPLE ALIGNMENTS FOR CONSERVED RNA STRUCTURES\n
Reads a set of maf files, calculates stats, scans with SISSIz or RNAz, outputs bed coordonates of high-confidence predictions/n
                    Usage:     java  ECSFinder.jar [options] -o output/directory -i input.maf (last parameter must be -i)/n
                    Output: 	Two types of results are produced:/n
                               (1)  the multiple sequence alignments associated to significant predictions/n
                                    are saved to files in the folder specified with the -o option./n
                                    File names correspond to their genomic coordinates in a .bed-compatible format. /n
                                    Ex: output/directory/chrX_12345_12500_80:75:23:14:8:z_300_+.aln/n
                               (2)  The genomic coordinates (.bed format) of ECSs are also written to the SDOUT/n
                                    (see additional options below)./n
                    output: Number chromosome start loci  end loci number of species MPI standard deviation Normalized shannon entropy/n
                            GC % content Gap %content SISSIZ Z-score/n
                               ***  N.B. the score field corresponds to the SISSIz Z-score x-100/n
                    Options:/n
                      -all             write out bed entires for all sampled alignments to STDOUT/n
                      -c     int       number of CPUs for calculations (default 4)/n
                      -g     int       max gap percentage of sequences for 2D prediction (default 50)/n
                      -sszr  double    report SISSIz+RIBOSUM hits below this Z-score (default -3)/n
                      -v               verbose (messy but detailed) output/n
                      -r                R-scape to be used/n

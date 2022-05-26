 **ECSFinder

**SCAN MULTIPLE ALIGNMENTS FOR CONSERVED RNA STRUCTURES <br />
Reads a set of maf files, calculates stats, scans with SISSIz or RNAz, outputs bed coordonates of high-confidence predictions<br />
                    Usage:     java  ECSFinder.jar [options] -o output/directory -i input.maf (last parameter must be -i)<br />
                    Output: 	Two types of results are produced:<br />
                               (1)  the multiple sequence alignments associated to significant predictions<br />
                                    are saved to files in the folder specified with the -o option.<br />
                                    File names correspond to their genomic coordinates in a .bed-compatible format. <br />
                                    Ex: output/directory/chrX_12345_12500_80:75:23:14:8:z_300_+.aln/n
                               (2)  The genomic coordinates (.bed format) of ECSs are also written to the SDOUT<br />
                                    (see additional options below).<br />
                    output: Number chromosome start loci  end loci number of species MPI standard deviation Normalized shannon entropy<br />
                            GC % content Gap %content SISSIZ Z-score<br />
                               ***  N.B. the score field corresponds to the SISSIz Z-score x-100<br />
                    Options:<br />
                      -all             write out bed entires for all sampled alignments to STDOUT<br />
                      -c     int       number of CPUs for calculations (default 4)<br />
                      -g     int       max gap percentage of sequences for 2D prediction (default 50)<br />
                      -sszr  double    report SISSIz+RIBOSUM hits below this Z-score (default -3)<br />
                      -v               verbose (messy but detailed) output<br />
                      -r                R-scape to be used<br />

 **ECSFinder

**SCAN MULTIPLE ALIGNMENTS FOR CONSERVED RNA STRUCTURES <br />
Reads a set of maf files, calculates stats, scans with SISSIz or RNAz, outputs bed coordonates of high-confidence predictions<br />
                    &ensp; Usage:     java  ECSFinder.jar [options] -o output/directory -i input.maf (last parameter must be -i)<br />
                    &ensp; Output: 	Two types of results are produced:<br />
                            &ensp;&ensp;   (1)  the multiple sequence alignments associated to significant predictions<br />
                            &ensp;&ensp;        are saved to files in the folder specified with the -o option.<br />
                            &ensp;&ensp;        File names correspond to their genomic coordinates in a .bed-compatible format. <br />
                            &ensp;&ensp;        Ex: output/directory/chrX_12345_12500_80:75:23:14:8:z_300_+.aln<br />
                            &ensp;&ensp;   (2)  The genomic coordinates (.bed format) of ECSs are also written to the SDOUT<br />
                            &ensp;&ensp;        (see additional options below).<br />
               &ensp;     output: Number chromosome start loci  end loci number of species MPI standard deviation Normalized shannon entropy<br />
                        &ensp;&ensp;    GC % content Gap %content SISSIZ Z-score<br />
                      &ensp;&ensp;         ***  N.B. the score field corresponds to the SISSIz Z-score x-100<br />
               &ensp;     Options:<br />
               &ensp;&ensp;       -all             write out bed entires for all sampled alignments to STDOUT<br />
               &ensp;&ensp;       -c     int       number of CPUs for calculations (default 4)<br />
               &ensp;&ensp;       -g     int       max gap percentage of sequences for 2D prediction (default 50)<br />
               &ensp;&ensp;       -sszr  double    report SISSIz+RIBOSUM hits below this Z-score (default -3)<br />
               &ensp;&ensp;       -v               verbose (messy but detailed) output<br />
               &ensp;&ensp;       -r                R-scape to be used<br />

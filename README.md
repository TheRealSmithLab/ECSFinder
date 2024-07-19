


<p align="center" width="100%">
    <img width="25%" src="https://user-images.githubusercontent.com/44384386/195381940-680064be-d53a-45b6-a5e1-a80ff1cb804e.jpg"> 
</p>

# ECSFinder

Scans multiple alignments for conserved RNA structures. Reads a set of maf files, calculates stats, scans with SISSIz, outputs bed coordonates of high-confidence predictions. We used the locally stable consensus secondary structure prediction algorithm [RNALalifold](https://www.tbi.univie.ac.at/RNA/RNALalifold.1.html) in a first pass to refine alignment boundaries before using [SISSIz](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-248) to assess if a predicted conserved structure is more likely than chance give the underlying alignments. The filtered alignments are then used in RNAalifold to get the minimal free energy of the structure and the pseudo-energy and [R-scape]([https://www.tbi.univie.ac.at/RNA/RNALalifold.1.html](http://eddylab.org/R-scape/)) is then used to calculate the minimal E-value of the helices and the number of significant base pairs. 


## Table of content

- [Installation](#installation)
    - [SISSIz](#sissiz)
    - [RNALalifold](#rnalalifold)
    - [ECSFinder](#ecsfinder)
    - [R-scape](#rscape)

- [Usage](#usage)
- [Output](#output)
- [Example](#example)
- [Links](#links)
    - Alignments
    - UCSC track hub

## Installation
Installation of SISSIz 2.0 is required for running ECSFinder. Version 1.0 can be found [here](https://github.com/ViennaRNA/SISSIz).

### SISSIz

To install SISSIz version 2.0:
```
cd /SISSIz
autoreconf -fvi
./configure
make
sudo make install (as root)
```
If you have no root rights on your system or prefer to install SISSIz
into a self-contained directory run configure for example like this:
```
./configure --prefix=/opt/programs/SISSIz --datadir=/opt/programs/SISSIz/share
```

### RNALalifold
Download the package on the ViennaRNA package [website](https://www.tbi.univie.ac.at/RNA/) and follow the [instructions](https://www.tbi.univie.ac.at/RNA/documentation.html#install)
```
tar -zxvf ViennaRNA-2.5.1.tar.gz
cd ViennaRNA-2.5.1
./configure
make
sudo make install
```

### ECSFinder
```
cd ECSFinder/src
javac ECSFinder.java
```
### R-scape
```
Download the source code [website]([https://www.tbi.univie.ac.at/RNA/](http://eddylab.org/R-scape/))
```
### ECSFinder
```
java ECSFinder [options] -o output/directory -i input.maf (last parameter must be -i, absolute path required)
 Options:
   -c int number of CPUs for calculations (default 4)
   -g int max gap percentage of sequences for 2D prediction (default 50)
   -sszr double report SISSIz+RIBOSUM hits below this Z-score (default -3)
   -v verbose (messy but detailed) output
```

## Output
Two outputs are produced: 
* A filtered MAF ready to be used in our pipeline (-output.maf)
* A text file containing all filtered out content (-removedLines.txt)

### ECSFinder
 Two types of results are produced:
   *  The multiple sequence alignments associated to significant predictions are saved to files in the folder specified with the -o option
      File names correspond to their genomic coordinates in a .bed-compatible format
      
      Ex: output/directory/X_9958021_9958096_11_92.2_0.1_0.16_66.8_0.4_304_+.aln
     
 ### File Name:
***
     
     1. Name of chromosome (X)
     2. Start loci (9958021)
     3. End loci (9958096)
     4. Number species (11)
     5. Mean pairwise identity (92.2)
     6. Standard deviation (0.1)
     7. Average Shannon entropy (0.16)
     8. GC content (66.8)
     9. Gap content (0.4)
    10. Z-score multiplied by -100 (304)
    11. Direction of the strand (+)
 ***   
     
   *  The genomic coordinates (.bed format) of ECSs are also written to the SDOUT
                                    
## Example
 ```

java ECSFinder -o output -c 10 -sszr -3.5 -i /home/vanda/46_mammals.epo.1 -output.maf

X	17713087	17713208	17:78.2:0.1:0.46:39.3:18.4	428	+
X	17713085	17713241	17:83.0:0.1:0.38:38.8:42.2	367	+
X	17713023	17713229	17:80.8:0.2:0.41:37.9:12.4	372	+
...
```

## Links

### Alignments
Alignments generated in our evolutionarily conserved structures scan can be downloaded [here](https://de.cyverse.org/data/ds/iplant/home/vandalovejoy/Data?type=folder&resourceId=84919f88-4984-11ed-9d98-90e2ba675364)

### UCSC track hub
To visualize the evolutionarily conserved structures using our UCSC track hub, you can enter this link: https://data.cyverse.org/dav-anon/iplant/home/vandalovejoy/ECS_hub.txt [here](https://genome.ucsc.edu/cgi-bin/hgHubConnect)
- ECS1.0 track contains the hits generated from our previous [study](https://academic.oup.com/nar/article/41/17/8220/2411364) including the trimmed consensus secondary structure for each alignmnent block.
- ECS2.0 track contains our new hits with the consensus secondary structure. The colors indicate that the ECS has one or more homologs of the same color.
- Homologs track contains the homologs of a subset of ECS2.0. Color is the same as that of it's ECS.

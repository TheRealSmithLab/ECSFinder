


<p align="center" width="100%">
    <img width="25%" src="https://user-images.githubusercontent.com/44384386/195381940-680064be-d53a-45b6-a5e1-a80ff1cb804e.jpg"> 
</p>

# ECSFinder

Scans multiple alignments for conserved RNA structures. Reads a set of maf files, calculates stats, scans with SISSIz or RNAz, outputs bed coordonates of high-confidence predictions. We used the locally stable consensus secondary structure prediction algorithm [RNALalifold](https://www.tbi.univie.ac.at/RNA/RNALalifold.1.html) in a first pass to refine alignment boundaries before using [SISSIz](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-248) to assess if a predicted conserved structure is more likely than chance give the underlying alignments.


## Table of content

- [Installation](#installation)
    - [SISSIz](#sissiz)
    - [RNALalifold](#rnalalifold)
    - [ECSFinder](#ecsfinder)
    - [MergeNFilter](#mergenfilter)
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
### MergeNFilter
It is essential to use MergeNFilter prior to using ECSFinder as it will prepare your MAF file for our pipeline. It will remove duplicate sequence, duplicate sequence IDs as well as ancestor sequences from it's input file.
```
javac MergeNFilter.java
```
## Usage
the input file needs to be in the same folder as MergeNFilter 
### MergeNFilter
```
java MergeNFilter 46_mammals.epo.1_19.maf
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
### MergeNFilter
Two outputs are produced: 
* A filtered MAF ready to be used in our pipeline
* A text file containing all filtered out content (_segmtl_dups.txt)

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
java MergeNFilter 46_mammals.epo.1_19.maf

java ECSFinder -o output -c 10 -sszr -3.5 -i /home/vanda/46_mammals.epo.1.maf

Stockholm file: alifold_0001.stk
Stockholm file: alifold_0002.stk
1       196390149       196390200       29      82.6    0.1     0.38    22.2    7.2     353     +
1       196386763       196386824       29      83.9    0.2     0.34    36.5    0.5     388     +
1       196386763       196386822       29      84.3    0.1     0.35    37.2    0.5     430     +
1       196386763       196386827       29      84.7    0.1     0.33    35.8    0.4     355     +
...
```

## Links

### Alignments
Alignments generated in our evolutionarily conserved structures scan can be downloaded [here](https://de.cyverse.org/data/ds/iplant/home/vandalovejoy/Data?type=folder&resourceId=84919f88-4984-11ed-9d98-90e2ba675364)

### UCSC track hub
To visualize the evolutionarily conserved structures using our UCSC track hub, click [here](https://genome.ucsc.edu/s/VandaLovejoy2/Paper%20ECS)
- ECS1.0 track contains the hits generated from our previous [study](https://academic.oup.com/nar/article/41/17/8220/2411364) including the trimmed consensus secondary structure for each alignmnent block.
- ECS2.0 track contains our new hits with the consensus secondary structure. The colors indicate that the ECS has one or more homologs of the same color.
- Homologs track contains the homologs of a subset of ECS2.0. Color is the same as that of it's ECS.

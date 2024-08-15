


<p align="center" width="100%">
    <img width="25%" src="https://user-images.githubusercontent.com/44384386/195381940-680064be-d53a-45b6-a5e1-a80ff1cb804e.jpg"> 
</p>

ECSFinder is a tool designed to scan multiple alignments for conserved RNA structures. It processes a set of MAF files, calculates key statistics, scans with SISSIz, and outputs BED coordinates of high-confidence predictions. The process begins by refining alignment boundaries using RNALalifold, which identifies locally stable RNA secondary structures. After this refinement, the alignments are analyzed with SISSIz to assess whether a predicted conserved structure is statistically more likely than expected by chance.

Following this, RNAalifold calculates the minimal free energy and pseudo-energy of the predicted structures, providing insights into their stability. R-scape is then used to evaluate the statistical significance of helices within the RNA structures, identifying significant base pairs that further validate the findings.

To enhance prediction accuracy, ECSFinder integrates a Generalized Linear Model (GLM) that predicts the likelihood of the identified RNA structures being true positives or false positives. This model considers features such as E-value, the number of significant base pairs, minimal free energy, pseudo-energy, and sequence conservation (MPI).

The result is a robust framework that not only identifies but also validates conserved RNA structures across multiple sequence alignments, providing output that can be visualized and further analyzed using genome browsers and other bioinformatics tools.


## Table of Contents

- [Installation](#installation)
  - [SISSIz](#sissiz)
  - [RNALalifold](#rnalalifold)
  - [ECSFinder](#ecsfinder)
  - [R-scape](#r-scape)
- [Usage](#usage)
- [Output](#output)
- [Example](#example)


## Installation
Installation of SISSIz 2.0 is required for running ECSFinder. Version 1.0 can be found [here](https://github.com/ViennaRNA/SISSIz).

### SISSIz

To install SISSIz version 2.0:
```sh
cd /SISSIz
autoreconf -fvi
./configure
make
sudo make install
```
If you have no root rights or prefer to install SISSIz into a self-contained directory, run configure like this:
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

Download the source code [website](http://eddylab.org/R-scape/)

## Usage
### ECSFinder
```
java ECSFinder [options] -o output/directory -i input.maf (last parameter must be -i, absolute path required)
 Options:
   -c int number of CPUs for calculations (default 4)
   -g int max gap percentage of sequences for 2D prediction (default 50)
   -sszr double report SISSIzhits below this Z-score (default -3)
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

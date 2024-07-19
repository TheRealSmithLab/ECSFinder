


<p align="center" width="100%">
    <img width="25%" src="https://user-images.githubusercontent.com/44384386/195381940-680064be-d53a-45b6-a5e1-a80ff1cb804e.jpg"> 
</p>

# ECSFinder

ECSFinder scans multiple alignments for conserved RNA structures. It processes a set of MAF files, calculates various statistics, scans with SISSIz, and outputs BED coordinates of high-confidence predictions. The tool integrates several key steps to ensure accurate identification and characterization of conserved RNA structures.

1. **Refinement of Alignment Boundaries**:
   - We use the locally stable consensus secondary structure prediction algorithm [RNALalifold](https://www.tbi.univie.ac.at/RNA/RNALalifold.1.html) in an initial pass to refine the boundaries of the alignments. This steps identified locally stable RNA secondary structures.

2. **Assessment of Conserved Structures**:
   - Next, the refined alignments are processed using [SISSIz](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-248). SISSIz evaluates whether a predicted conserved structure is statistically more likely than expected by chance, given the underlying sequence alignments. This statistical assessment using the Z-score, helps in identifying potentially significant RNA structures.

3. **Energy Calculations**:
   - The filtered alignments are then analyzed with RNAalifold to determine the minimal free energy of the predicted structures. RNAalifold also calculates the pseudo-energy, which provides additional insights into the stability and likelihood of the RNA structures.

4. **Helix Significance Calculation**:
   - Finally, [R-scape](http://eddylab.org/R-scape/) is employed to calculate the minimal E-value of the helices within the RNA structures and the number of significant base pairs. R-scape's analysis helps in identifying helices that are statistically significant, further validating the conserved RNA structures.

5. **Prediction Model with Generalized Linear Model (GLM)**:
   - ECSFinder integrates a Generalized Linear Model (GLM) in `prediction_GLM.R` to predict the likelihood of the identified RNA structures being true positives (TP) or false positives (FP). The model uses several features extracted during the analysis, including:
     - `log10(E-value)`: Logarithm of the minimal E-value from the helices from R-scape.
     - `number of significant base pairs`: Number of covarying base pairs from R-scape.
     - `MFE`: Minimal free energy from RNAalifold.
     - `pseudo_energy`: Pseudo-energy from RNAalifold.
     - `MPI` (Mean Pairwise Identity): A measure of sequence conservation across the alignment.
   - The GLM uses these features to provide a probabilistic prediction, aiding researchers in distinguishing between likely true and false conserved RNA structures.

By combining these tools and steps, ECSFinder provides a robust framework for the identification and analysis of conserved RNA structures across multiple sequence alignments. The tool generates BED coordinates for high-confidence predictions, which can be visualized and further analyzed using genome browsers and other bioinformatics tools. The integration of the GLM enhances the accuracy and reliability of the predictions, providing an additional layer of validation for the identified RNA structures.


## Table of content

## Table of Contents

- [Installation](#installation)
  - [SISSIz](#sissiz)
  - [RNALalifold](#rnalalifold)
  - [ECSFinder](#ecsfinder)
  - [R-scape](#rscape)
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
```
Download the source code [website]([https://www.tbi.univie.ac.at/RNA/](http://eddylab.org/R-scape/))
```
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

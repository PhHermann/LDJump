# LDJump
LD Jump is an R Package to estimate variable recombination rates from population genetic data. 
It is a unix based program (with a necessary installation of LDHat (see <https://github.com/auton1/LDhat>)), able to estimate the recombination map of sequences in fasta and vcf format. 
First, the sequences are divided in short segments of user defined length. The recombination rate is estimated for every segment with a regression model. 
This set of estimates is fed in a segmentation algorithm (SMUCE) to estimate the breakpoints of the recombination landscape.

### Author (Requests)
Please contact me in case of questions, comments, bug reports, etc... 
Author: Philipp Hermann
E-Mail: <philipp.hermann@jku.at>

## Dependencies
This unix package makes use of several functions of other R-packages listed as follows: 

* Unix Operating System
* R (>= 2.10)
* adegenet (>= 2.0.1)
* ape
* genetics (>= 1.3.8.1)
* Biostrings (>= 2.38.4)
* stepR
* seqinr (>= 3.1-3)
* quantreg
* vcfR (>= 1.5.0)

The uploaded version on CRAN will install these packages automatically. 

## Installation
The ZIP-File of the package can be downloaded via "Clone or download". The installation is performed as usually in R via a the command: 

    install.packages("/PathToLDJump/LDJump_<version>.tar.gz", repos=NULL, type="source")

## Usage

After loading the Package in the Workspace one can use the main function *LDJump* in order to estimate the variable recombination rate for the population of interest. We recommend to use **LDJump** with the population in *fasta*-Format. Alternatively, restructuring of *vcf*-Files is also possible with a implemented function *vcfR_to_fasta*. : 

    require(LDJump)
    LDJump(seqName, alpha = 0.05, segLength = 1000, pathLDHat = "", format = "fasta", refName = NULL, start = NULL)
    
Detailed descriptions of the main functions and all adjacent functions computing the recombination map can be found via 

    ?LDJump
    
An example with files can be found in **Examples**. 

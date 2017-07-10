# LDJump
LD Jump is an R Package to estimate variable recombination rates from population genetic data. 
It is a unix based program (with a necessary installation of LDHat (see https://github.com/auton1/LDhat)), able to estimate the recombination map of sequences in fasta and vcf format. 
First, the sequences are divided in short segments of user defined length. The recombination rate is estimated for every segment with a regression model. 
This set of estimates is fed in a segmentation algorithm (SMUCE) to estimate the breakpoints of the recombination landscape.

## Dependencies

* item Test

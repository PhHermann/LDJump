# Example

Here, we upload all necessary files for the example of an application of **LDJump**. 

* TSI_Sample.zip: This file contains the sample of 50 invidivuals of the Toscan (Italian) population of the 1000 Genomes project. 
* Lookups.zip: This file contains the lookup table of 100 sequences with a theta of 0.005. 

The required R-command for the estimation of the recombination map with **LDJump** is the following: 

    require(LDJump)
    results = LDJump("/pathToSample/TSI_CH21_41187000_41290679.fa", alpha = 0.05, segLength = 1000, pathLDHat = "/pathToLDHat",
                format = "fasta", refName = NULL, thth = 0.005)
    plot(results)
    
With the *plot*-Function of the package *stepR* one will obtain the estimated map with the estimated recombination rates plotted on the y-axis and the according segment number on the x-axis. 

Note: Both ZIP-Files have to be unzipped before usage in the R-function *LDJump*. 


# Example

Here, we upload all necessary files for the example of an application of **LDJump**. 

* [SimulatedPopulations.zip](https://github.com/PhHermann/LDJump/blob/master/Example/SimulatedPopulations.zip): This file contains three different simulation setups of 16 individuals (=sequences) with background rates of 0.001, 0.0054, and 0.001. This setup contains 15 *Hotspots* of 8 to 40 folds of the background rates evenly distributed across the sequence. 
<!-- TSI_Sample.zip This file contains the sample of 50 invidivuals of the Toscan (Italian) population of the 1000 Genomes project-->

The required R-command for the estimation of the recombination map with **LDJump** for the recombination map with a background rate of 0.054 is the following: 

```R
require(LDJump)
results = LDJump("/pathToSample/HatLandscapeN16Len1000000Nrhs15_th0.01_540_1.fa", alpha = 0.05, segLength = 1000, 
pathLDhat = "/pathToLDhat", pathPhi = "/pathToPhi", format = "fasta", refName = NULL)
postscript("Results.eps", horiz = F)
plot(results[[1]], xlab = "Segments", ylab = "Estimated Recombination Rate", main = "Estimated recombination map with LDJump")
dev.off()

```

    
With the *plot*-Function of the package *stepR* one will obtain the estimated map with the estimated recombination rates plotted on the y-axis and the according segment number on the x-axis. The *postscript* and *dev.off* commands before and after plotting the results, respectively, will save the result in EPS-format. 

**LDJump** returns a list of 7 elements which contains the estimated recombination map, the constant recombination rate estimates per segment, the calculated summary statistics (in a matrix), the type I error probability, the sample size, the sequence length and the segment lenght used. 
For a constant recombination rate estimation only the latter six elements are returned. 

Notes: 
* Both ZIP-Files have to be unzipped before usage in the R-function *LDJump*. 
* Runtime of one of these examples will be at least one hour. 
* We recommend when *LDJump* is applied with different populations/examples simultaneously to create directories and start *LDJump* separately from these directories. 


# Example

Here, we upload all necessary files for the example of an application of **LDJump**. 

* [SimulatedPopulations.zip](https://github.com/PhHermann/LDJump/blob/master/Example/SimulatedPopulations.zip): This file contains three different simulation setups of 16 individuals (=sequences) with background rates of 0.001, 0.0054, and 0.001. This setup contains 15 *Hotspots* of 8 to 40 folds of the background rates evenly distributed across the sequence. 
<TSI_Sample.zip This file contains the sample of 50 invidivuals of the Toscan (Italian) population of the 1000 Genomes project>
* [Lookups.zip](https://github.com/PhHermann/LDJump/blob/master/Example/Lookups.zip): This file contains the lookup table of 100 sequences with a theta of 0.005. 

The required R-command for the estimation of the recombination map with **LDJump** for the recombination map with a background rate of 0.054 is the following: 

```R
require(LDJump)
results = LDJump("/pathToSample/HatLandscapeN16Len1000000Nrhs15_th0.01_540_1.fa", alpha = 0.05, segLength = 1000, 
pathLDHat = "/pathToLDHat", format = "fasta", refName = NULL, thth = 0.005)
plot(results)
```
    
With the *plot*-Function of the package *stepR* one will obtain the estimated map with the estimated recombination rates plotted on the y-axis and the according segment number on the x-axis. 

Note: Both ZIP-Files have to be unzipped before usage in the R-function *LDJump*. 


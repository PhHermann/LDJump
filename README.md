# LDJump
**LDJump** is an R package to estimate variable recombination rates from population genetic data. 
It is a unix based program (with a necessary installation of [PhiPack](<https://www.maths.otago.ac.nz/~dbryant/software.html>)), able to estimate the recombination map of sequences in fasta and vcf format. 
First, the sequences are divided in short segments of user defined length. The constant recombination rate is estimated for every segment with a regression model. 
This set of estimates is fed in a segmentation algorithm (SMUCE) to estimate the breakpoints of the recombination landscape. A [PDF Manual](./LDJump.pdf) with complete documentation of each function is also available. The publication introducing **LDJump** can be found [here](<https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12994?af=R>).

### Author (Requests)
Please contact me in case of questions, comments, bug reports, etc...

    Author: Philipp Hermann
    E-Mail: philipp.hermann@jku.at

## Dependencies & System Requirements
This package makes use of several functions of other R-packages, of the software package [PhiPack](<https://www.maths.otago.ac.nz/~dbryant/software.html>), and it can optionally also make use of [LDhat](<https://github.com/auton1/LDhat>) in order to decrease the runtime of **LDJump**. Here we provide a list of used R-packages, as well as, other software packages: 

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
* PhiPack
* LDhat (2.2, optional -> to decrease runtime)
* snow
* scrm (1.7.1, optional to train the regression model under different scenarios than considered)
* ms2dna (1.16, optional to train the regression model under different scenarios than considered)

Notice that **Biostrings** needs to be installed via [Bioconductor](<http://bioconductor.org/packages/release/bioc/html/Biostrings.html>).  [PhiPack](<https://www.maths.otago.ac.nz/~dbryant/software.html>) as well as the functions *dos2unix* and *awk* neeed to be installed too. Additionally, we recommend to also install [LDhat](<https://github.com/auton1/LDhat>) which will then enable to compute estimates at lower computational cost. The software packages [scrm](<https://github.com/scrm/scrm>) and [ms2dna](<http://guanine.evolbio.mpg.de/bioBox/>)) only need to be installed when the user wants to estimate a regression model based on simulations under a user input demographic scenario, see [Update to Version 0.2.2](<###Update to Version 0.2.2>).


## Installation
The ZIP-File of the package can be downloaded via "Clone or download" as well as the command line: 

```markdown
git clone https://github.com/PhHermann/LDJump.git
R CMD build LDJump
R CMD INSTALL LDJump_<version>.tar.gz
``` 

Using *R* one can also install the *.tar.gz* file with the following command: 
```R
install.packages("/PathToLDJump/LDJump_<version>.tar.gz", repos=NULL, type="source")
```

## Usage

After loading the package in the workspace one can use the main function *LDJump* in order to estimate the variable recombination rate for the population of interest. We recommend to use **LDJump** with the population in *fasta*-Format. Alternatively, restructuring of *vcf*-Files is also possible with a implemented function *vcfR_to_fasta*. Therefore, we used the reference sequence from <http://phase3browser.1000genomes.org/Homo_sapiens> for our example: 

```R
require(LDJump)
LDJump(seqName, alpha = 0.05, segLength = 1000, pathLDhat = "", pathPhi = "", format = "fasta", refName = NULL, 
      start = NULL, constant = F, status = T, cores = 1, accept = F, demography = F, out = "")
```

Detailed descriptions of the main functions and all adjacent functions computing the recombination map can be found via e.g.

```R
?LDJump
```

We provide examples with files in *[Example](./Example)*. Previous versions (before version 0.2.1) required a set of lookup-tables of [LDhat](<https://github.com/auton1/LDhat>) which can still be found in *[Lookup Tables](./Lookups)*. However, we recommend to update to the most recent version of **LDJump**. 

A full path to the *Phi* file of [PhiPack](<https://www.maths.otago.ac.nz/~dbryant/software.html>) needs to be provided as follows *pathPhi = "/path/to/Phi"*. In order to use [LDhat](<https://github.com/auton1/LDhat>) to quickly calculate some of the summary statistics please set *pathLDhat = "/path/to/LDhat-master"*.

**LDJump** can also be used under a set of type I error probabilities alpha. Therefore, the parameter *alpha* can also be fed with a vector of values such as:

```R
require(LDJump)
LDJump(seqName, alpha = c(0.1, 0.05, 0.01), segLength = 1000, pathLDhat = "", pathPhi = "", format = "fasta",
       refName = NULL, start = NULL, constant = F, status = T, cores = 1)
```

**LDJump** is designed to estimate recombination rates from segments containing information (by SNPs). Therefore, the program checks all segments (based on the given segment lengths) for the number of SNPS. In case of segments without SNPs, the program will inform the user and ask for input providing the following two options: 
* "n": **LDJump** will be stopped and the user should retry with a larger segment length. 
* "y": **LDJump** will continue and impute the recombination rate of the segment without SNPs with a weighted mean of adjacent segments. 

In order to avoid the user input and accept that **LDJump** imputes recombination rates for these segments, we have added the parameter *accept* (by default *FALSE*), which should then be set to *TRUE*. 

We also included a logical parameter *constant* in **LDJump**, which is *FALSE* by default to estimate **variable recombination rates**. In the case that *constant* is set to *TRUE*, **LDJump** will provide a **constant recombination rate estimator** of the whole sequence. 

A logical parameter *rescale* enables to transform the sequence positions to the unit interval if set to *TRUE*.

An integer parameter *cores* enables to parallelize **LDJump**, where the parameter equals the number of cores on which **LDJump** should run. 

A logical parameter *status*, which is *TRUE* by default, prints the current status of the calculated segment on screen or to the file *LDJump_Status.txt* in case of applying **LDJump** on several cores in parallel. 

A logical parameter *demography*, which is *FALSE* by default, enables to estimate the recombination rates using a generalized additive regression model (GAM) which is based on simulated samples of populations undergoing demographic effects. 

### Update to Version 0.2.2
We have implemented a function *calcRegMod* which computes the regression model for constant estimates by simulating under user defined populations. Therefore, we provide the following example how to obtain the regression  models for populations under a bottleneck followed by rapid growth in population sizes. 

```R
require(LDJump)
scenario =  " -eG 0.0 0 -eG 0.42 -100 -eG 0.5 100 "
simulatedData = calcRegMod(nsim=100,pathToScrm=pathToScrm,scenario=scenario,pathToMs2dna=pathToMs2dna, status = T, pathLDhat = "/path/To/LDhat", pathPhi = "/path/To/Phi")
regMod = simulatedData[[1]]
result = LDJump(seqName = "path/To/seqName", pathLDhat = "/path/To/LDhat", pathPhi = "/path/To/Phi", segLength = 1000, alpha = 0.05, status = T, demography = F, regMod = regMod, status = T, accept = T, cores = 1)
```

Notice that any demographic scenario which is in the range of [scrm](<https://github.com/scrm/scrm/wiki/Command-Line-Options>) can be provided with the simulation function. **However, we want to stress that the resulting model is not checked for its underlying assumptions of the residuals. Nevertheless, this can and should be done manually, given that the model is in hand after the simulations are performed.** Moreover, one can also adapt the sample sizes of the populations (n), the sequence lengths of the simulated populations (len) and the recombination rates. The default values are listed as follows. 
```R
n = c(10, 16, 20)
len = c(500, 1000, 2000, 3000, 5000)
```

For each combination of sample sizes and sequence lengths 100 recombination rates are simulated (uniformly distributed per interval) for the 5 following intervals: [0,0.01], [0.01,0.02], [0.02, 0.05], [0.05, 0.1], and [0.1, 0.2]. As noted above, we added the option to provide a vector of recombination rates under which the recombination map should be simulated. Then the mentioned setup using intervals is ignored. 

Although we did not explore **LDJump** under other scenarios such as population structure, we note that principially (given its implementation in [scrm](<https://github.com/scrm/scrm/wiki/Command-Line-Options>) it is possible to simulate under any scenario using our *calcRegMod* function. **However, one has to be careful that we did not address other scenarios than the demographic example as provided with the R-package. Hence, one has to be careful with LDJump and self-simulated regression models.**

We added an optional parameter *out* which is a prefix for all outfiles created within **LDJump**. By default it is an empty string, but when changed to any user defined-string, one can run **LDJump** in parallel on several data sets from the same working directory without interfering output files. 

### Recommendations
We recommend to run **LDJump** from the same path where the sample file is located in order that all created temp files will be deleted after completion. 

# Reference 
Hermann, P., Heissl, A., Tiemann‐Boege, I., and Futschik, A. (2019), LDJump: Estimating Variable Recombination Rates from Population Genetic Data. Mol Ecol Resour. [doi:10.1111/1755-0998.12994](<https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12994?af=R>). 

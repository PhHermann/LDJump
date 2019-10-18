# LDJump
**LDJump** is an R package to estimate variable recombination rates from population genetic data. 
It is a unix based program (with a necessary installation of [PhiPack](<https://www.maths.otago.ac.nz/~dbryant/software.html>)), able to estimate the recombination map of sequences in fasta and vcf format. 
First, the sequences are divided in short segments of user defined length. The constant recombination rate is estimated for every segment with a regression model. 
This set of estimates is fed in a segmentation algorithm (SMUCE) to estimate the breakpoints of the recombination landscape. A [PDF Manual](./LDJump.pdf) with complete documentation of each function is also available. The publication introducing **LDJump** can be found [here](<https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12994?af=R>).

### Author (Requests)
Please contact me in case of questions, comments, bug reports, etc...

    Author: Philipp Hermann, Fardokhtsadat Mohammadi
    E-Mail: philipp.hermann@jku.at, fardokht.fm@gmail.com

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
* parallel
* [VCFTools](<https://vcftools.github.io/index.html>) (only required for VCF-files)

Notice that **Biostrings** needs to be installed via [Bioconductor](<http://bioconductor.org/packages/release/bioc/html/Biostrings.html>).  [PhiPack](<https://www.maths.otago.ac.nz/~dbryant/software.html>) as well as the functions *dos2unix* and *awk* neeed to be installed too. Additionally, we recommend to also install [LDhat](<https://github.com/auton1/LDhat>) which will then enable to compute estimates at lower computational cost. The software packages [scrm](<https://github.com/scrm/scrm>) and [ms2dna](<http://guanine.evolbio.mpg.de/bioBox/>)) only need to be installed when the user wants to estimate a regression model based on simulations under a user input demographic scenario, see [Update to Version 0.2.2](<###Update to Version 0.2.2>). In regard to [Update to Version 0.3.1](<###Update to Version 0.3.1>), **VCFTools** is a necessary package. For directions on how to install it, please refer to their [package documentation](<https://vcftools.github.io/examples.html>).


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

### Update to Version 0.3.1
We have extended LDJump in this update to handle VCF-files successfully. The most important point to adress in regard to working with VCF-files is that converting huge VCF-files to FASTA at once using g.e. functions like *vcfR2DNAbin* from **vcfR** proves to be very difficult and time consuming. The workflow of LDJump, nevertheless, requires FASTA files to compute the desired summary statistics, which is why we decided to first slice the VCF-file into several segments using [VCFTools](<https://vcftools.github.io/index.html>). The segmentation is the very first step when running **LDJump** on VCF-files. 

In order to run LDJump on VCF-files, two types of files are required:
* VCF-file: The VCF-file which will be used for the analysis.
* Reference-file: A fasta file of the very same sequence range used to convert the VCF-file to FASTA. 

Additionally to the filename/filepath, the command to run **LDJump** on VCF-files also requires the chromosome number, starting and ending position. For example, we could be running our analysis on chromosome 21 - ranging from base pair 41187000 to 41290679. Possible R commands, then, would look like this: 

```R
setwd(~/your_directory) #consider setting your working directory in the folder 
                        #in which you have all the files that you will be using (VCF-file, Reference-File)
                         
vcf_file = "TSI_21_41187000_41290679.vcf"
ref_seq = "Reference_CH21_41187000_41290679.fasta"

startofseq = 41187000
endofseq = 41290679
chr = 21
```

A possible command on running LDJump looks like this: 

```R
LDJump(vcf_file, chr = chr , segLength = 1000, cores = 8, pathPhi = "/path/To/Phi", format = "vcf", refName = ref_seq, lengthofseq=10000, startofseq = startofseq, endofseq = endofseq)
```

As it is mentioned above, other parameters such as *segLength, cores, pathPhi* are being used, which are generally necessary and can be referred to in the usual [LDJump documentation](./LDJump.pdf). 

The workflow for FASTA remained the same, but for VCF-files it has changed:

1. In *LDJump.R*, it checks whether your file is a VCF-File. 
   If yes, it will create a temporary folder "temp" in which the converted files will be saved.
2. If format == "vcf", the program redirects into *vcf_statistics.R*. 
   Here, it will initialize some variables that will be needed for our main looping and slicing explained in step 3.
3. In the looping, we segment the vcf-file according to the segment-range we set in the parameters (g.e. segLength = 1000). 
   Then, we convert it using *vcfR_to_fasta.R*.
4. *vcfR_to_fasta.R* reads into the newly segmented VCF-file and converts it into DNABin. 
   This conversion is only possible because we are first creating a subset of our Reference-Fasta-Sequence. 
   Each subsetted Fasta-sequence will then be used to convert the current VCF-segment.
   Afterwards, we write out the newly created FASTA-file into our temporary folder. 
5. The next step runs according to the usual procedure in LDJump.
   We use sapply to iterate over all the segments and to compute our summary_statistics,
   which are then combined to a data frame called "helper".

Last but not least: Before LDJump is executed again, the "temp" folder has to be deleted in order for LDJump to work.


**For any further questions, please contact Fardokhtsadat Mohammadi: fardokht.fm@gmail.com**



### Recommendations
We recommend to run **LDJump** from the same path where the sample file is located in order that all created temp files will be deleted after completion. 

# Reference 
Hermann, P., Heissl, A., Tiemann‐Boege, I., and Futschik, A. (2019), LDJump: Estimating Variable Recombination Rates from Population Genetic Data. Mol Ecol Resour. [doi:10.1111/1755-0998.12994](<https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12994?af=R>). 

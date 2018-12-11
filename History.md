## Source History

The most recent version of LDJump is 0.2.2 For installation of **LDJump** please see [README.md](README.md) and can be downloaded in this folder. 

### Update Version 0.1.1.
The functions *LDJump* and *getSmuce* are extended with a parameter *constant*. 
By default this is set to *FALSE* such that variable recombination rates are estimated. 
If set to *TRUE*, then a constant recombination rate is estimated for the whole sequence. 

### Update Version 0.1.2
The *smuceR* function to apply SMUCE is exchanged with the function *stepFit* due to the update of the stepR package and the notification that the first function might be removed in future updates. 
Moreover, a note due to package imports and dependencies has been fixed. 

### Update Version 0.1.3
We have changed several function calls due to package imports and dependencies in the previous version. 
We also updated the main function in order to be able to efficiently apply *LDJump* under a vector of type I error probabilities. Moreover, the output has more detailed descriptions. 
We also added an optional parameter to rescale the sequence length to the unit interval. 

### Update Version 0.1.4
Fixed minor bugs and paths to LDhat. Also an optional parameter is added which prints the current status of the calculated segment on screen. Example files are also updated. 

### Update Version 0.1.5
Added the option to impute recombination rates of segments having less than 2 SNPs. User input is required to continue the computation with **LDJump**. 

### Update Version 0.1.6
Fixed minor bug and added option to *fgt_rrate_dpr.R*. 

### Update Version 0.1.7
**LDJump** can now be applied in parallel on several cores. Therefore, an additional parameter is added to the command which reflects the integer number of cores on which the recombination map should be estimated. 

### Update Version 0.1.8
Thanks to user request we have added the option to ignore the user input in case of segments lacking SNPs. Moreover, slight adaptions to commands using the working directory were made. 

### Update Version 0.1.9
We have also simulated under demographic effects and provide a regression model for this setup. This is implemented via an logical parameter *demography*, which will use this regression model if the recombination rates should be estimated under demographic effects. 

### Update Version 0.1.10
Fix issues with paths. *LDJump* should be applied with the full path of the file containing the sequences  (\code{seqName}). Moreover, the function *impute* was rewritten. 

### Update Version 0.2.1
We have implemented the function *calcRegMod* in order to estimate the regression function for constant estimates under simulated populations with a user input demographic scenario. This function requires the two software packages [scrm](<https://github.com/scrm/scrm>) and [ms2dna](<http://guanine.evolbio.mpg.de/bioBox/ms2dna_1.16.tgz>) to be installed. Therefore, we also included an optional parameter in the main function in oder to apply **LDJump** under the newly trained regression model. 

We added the parameter *out* which enables to specify a string added to the output-files. Hence, it is not required anymore to run *LDJump* from different directories when applying it on several data sets in parallel. 

### Update Version 0.2.2
We updated the regression models for constant recombination rates in **LDJump** and therefore, also updated all other models such as the bias correction model or the demography model. These models are based on summary statistics computed from among others the software package [PhiPack](<https://www.maths.otago.ac.nz/~dbryant/software/>). The computing time of **LDJump** is reduced drastically. Additionally, one can still use [LDhat](<https://github.com/auton1/LDhat>) in order to reduce the runtime of **LDJump**. 


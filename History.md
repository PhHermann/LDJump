## Source History

The most recent version of LDJump is 0.1.8. For installation of **LDJump** please see [README.md](README.md) and can be downloaded in this folder. 

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

### Update Version 0.1.9
Fix issues with paths. *LDJump* should be applied with the full path of the file containing the sequences  (\code{seqName}). 


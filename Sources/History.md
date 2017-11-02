## Source History

The most recent version of LDJump is 0.1.3 and can be downloaded in this folder. 

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

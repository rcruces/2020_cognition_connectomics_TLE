# Code  

| **File**  | **Description**  |
|---|---|
| [`connectome_1_pre`](https://github.com/rcruces/2020_cognition_connectomics_TLE/blob/master/code/connectome_1_pre)  |  Prepares the necessary files to run ACT and SIFT on the already corrected DWI |
| [`connectome_2_sift`](https://github.com/rcruces/2020_cognition_connectomics_TLE/blob/master/code/connectome_2_sift)  | Runs ACT then SIFT  |
| [`connectome_3_custom_seg`](https://github.com/rcruces/2020_cognition_connectomics_TLE/blob/master/code/connectome_3_custom_seg)  | Prepares the Nodes NIFTI file by merging Dextrieux atlas and volbrain  |
| [`rCCA_cognition-connectomics.R`](https://github.com/rcruces/2020_cognition_connectomics_TLE/blob/master/code/rCCA_cognition-connectomics.R)  | Example of  regularized Canonical Correlations Analysis  |
| [`functions_cca.R`](https://github.com/rcruces/2020_cognition_connectomics_TLE/blob/master/code/functions_cca.R)  | Functions for the regularized Canonical Correlations Analysis  |
| [`functions_extra.R`]()  | Additional functions  |

The [`R`](https://www.r-project.org) code was test and run on version 3.4.4 and is actually being tested in the lastest version 3.6.2.  
Additionally some R packages are required to be installed.  

``` R
install.packages(CCA)
install.packages(scales)
install.packages(candisc)
install.packages(boot)
install.packages(foreach)
install.packages(doMC)
````

# Multidimensional Associations Between Cognition and Connectome Organization in Temporal Lobe Epilepsy
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/4293e4fe51114d02bb930c7d4b28f78b)](https://www.codacy.com/gh/rcruces/2020_cognition_connectomics_TLE/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=rcruces/2020_cognition_connectomics_TLE&amp;utm_campaign=Badge_Grade)  

Raúl Rodríguez-Cruces  [raul.rodriguezcruces@mcgill.ca](mailto:raul.rodriguezcruces@mcgill.ca) & [raulrcruces@gmail.com](mailto:raulrcruces@gmail.com)  
Boris C. Bernhardt  [boris.bernhardt@mcgill.ca](mailto:boris.bernhardt@mcgill.ca)  
Luis Concha [lconcha@unam.mx](mailto:lconcha@unam.mx)  

**Cite:** Rodríguez-Cruces, R., Bernhardt, B. C., & Concha, L. (2020). Multidimensional associations between cognition and connectome organization in temporal lobe epilepsy. NeuroImage, 116706.  

**DOI:** [j.neuroimage.2020.116706](https://doi.org/10.1016/j.neuroimage.2020.116706).  
**OSF:**  
**Preprint:** avaliable at [bioRxiv doi: https://doi.org/10.1101/675884](https://www.biorxiv.org/content/10.1101/675884v3.abstract).   
    
## Abstract
**Objective.** Temporal lobe epilepsy (TLE) is known to affect large-scale structural networks and cognitive function in multiple domains. The study of complex relations between structural network organization and cognition requires comprehensive analytical methods and a shift towards multivariate techniques. Here, we sought to identify multidimensional associations between cognitive performance and structural network topology in TLE.  
  
**Methods.** We studied 34 drug-resistant adult TLE patients and 25 age- and sex-matched healthy controls. Participants underwent a comprehensive neurocognitive battery and multimodal MRI, allowing for large-scale connectomics, and morphological evaluation of subcortical and neocortical regions. Using canonical correlation analysis, we identified a multivariate mode that links cognitive performance to a brain structural network. Our approach was complemented by bootstrap-based hierarchical clustering to derive cognitive subtypes and associated patterns of macroscale connectome anomalies.  
  
**Results.** Both methodologies provided converging evidence for a close coupling between cognitive impairments across multiple domains and large-scale structural network compromise. Cognitive classes presented with an increasing gradient of abnormalities (increasing cortical and subcortical atrophy and less efficient white matter connectome organization in patients with increasing degrees of cognitive impairments). Notably, network topology characterized cognitive performance better than morphometric measures did.  
  
**Conclusions.** Our multivariate approach emphasized a close coupling of cognitive dysfunction and large-scale network anomalies in TLE. Our findings contribute to understand the complexity of structural connectivity regulating the heterogeneous cognitive deficits found in epilepsy.  
  
 ## Repository content
 | Directories   | Description                                                                                                                                                                                                                                                                             |
|---------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [`./code`](https://github.com/rcruces/2020_cognition_connectomics_TLE/tree/master/code)      | `bash` functions for the connectome processing and R code for the rCCA and BASC analysis.                                                                                                                                                                                               |
| [`./databases`](https://github.com/rcruces/2020_cognition_connectomics_TLE/tree/master/databases) | Databases and tables necessary to run the code.                                                                                                                                                                                                                                         |
| [`./thickness`](https://github.com/rcruces/2020_cognition_connectomics_TLE/tree/master/thickness) | Contains all the surface files for each subject that was used in this study. The name is `ID_side_fsaverage5_20.mgh`, where ID is the identification number, side is the hemisphere (lh for left and rh for right), fsaverage5 is the surface used and 20 is the FWHM at 20 milimiters. |

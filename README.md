# RTGP
Code for Bayesian scalar-on-image regression with a relaxed-thresholded Gaussian process prior.

![alt text](https://github.com/annamenacher/RTGP/blob/main/surface_beta_map.png)


## Data

### Abstract
- Simulation study:
  - Volumetric-based representation: 2D images generated via data generating process described in Section C.3 of appendix in thesis.
  - Surface-based representation: Subsampled surface images generated via data generating process described in Section 4.3.1.1 in thesis.
- ABCD study: Surface-based representation of test statistic values from first-level fMRI analysis of emotional n-back task with 2 vs 0 back contrast.

- X: NxM matrix: vectorised images of test statistics of task fMRI
- X_confounds: NxP matrix: confounding variables
- y: Nx1 vector: intelligence scores

### Availability 
- Simulation study: Data can be generated via simulation study script '.../simulation_study/volume/data_circle.R' for the volumetric-based simulation study and '.../simulation_study/surface/data_generation.R' for the surface-based simulation study. We additionally provide an example simulation study containing 10 datasets (= 2D images with circle as effect of interest in the volume-based setting and sub-sampled surface maps in the surface-based setting) generated with the setting of sample size N=500 and SNR = 1/0.2 to showcase the reproducibility of our method. 
- ABCD study: Healthcare data requires access permission and cannot be attached here (https://abcdstudy.org/).

### Description
Simulated 2D volumetric images provided in this repositry under .../simulation_study/volume/N500_poly10_signal1_noise0.2_effect124/, simulated surface images provided in this repositry under .../simulation_study/surface/N500_L100_phi0.077_nu2/, and surface cifti files for COPE and standard deviation of COPE from first-level analysis provided from the ABCD study (not uploaded due to access restrictions).

## Code

### Abstract 
Code to run all baseline models, such as the frequentist approaches, Ridge regression and LASSO regression, the Bayesian regression models that take all the spatial locations as inputs, BR + Normal and BR + Horseshoe, and the Gaussian process regressions, GPR + Normal and GPR + Horseshoe, and our models RTGP-VI and RTGP-Gibbs. Additionally, we provide all code to generate figures & tables.

### Description 

Simualtion study (volume):
- **data_circle.R** : Generates simulation study data.
- **simulation_study_volume.R** : Performs parameter estimation and inference for models, such as RTGP-VI, RTGP-Gibbs, STGP, and GPR + Normal / Horeshoe.

Simulation study (surface):
- **data_generation.R** :
- **distance_matrix_left.py** :
- **distance_matrix_right.py** :
- **eigendecomposition_left_sphere.R** :
- **eigendecomposition_right_sphere.R** :
- **evaluate_datasets.R** :
- **evaluate_simulation_study_surface_ELBO.R** :
- **evaluate_simulation_study_surface_hyperparameters.R** :
- **evaluation_data_removal.R** :
- **simstudy_surface.R** :
- **simstudy_surface_data_removal.R** :
- **simstudy_surface_sensitivity_variance_hyperparameters.R** :
  
ABCD study:
- **BR_Normal_BR_Horseshoe.R** :
- **GPR_Normal_GPR_Horseshoe.R** :
- **RTGP_VI_left_hemisphere.R** :
- **RTGP_VI_right_hemisphere.R** :
- **RTGP_VI_left_hemisphere_CV.R** :
- **RTGP_VI_right_hemisphere_CV.R** :
- **Ridge_LASSO.R** :
- **analysis_eigendecomposition.R** :
- **brainsmash_distance_matrix_left_hemisphere.py** :
- **brainsmash_distance_matrix_right_hemisphere.py** :
- **data_preprocessing_step1.R** :
- **data_preprocessing_step2.R** :
- **eigendecomposition_left_cortex_L800_CV.R** :
- **eigendecomposition_right_cortex_L800_CV.R** :
- **evaluate_CV.R** :
- **evaluation_GPR.R** :
- **evaluation_RTGP.R** :
- **massunivariate_analysis.R** :

### Optional Information 

**R libraries:**

The following R packages are necessary to successfully use the codes:

- MASS
- BayesGPfit
- parallel
- fields
- ggplot2
- reshape2
- gridExtra
- truncnorm
- fastBayesReg
- STGP
- ggpointdensity
- xtable
- data.table
- ciftiTools
- bayestestR


Additional: 
- Software package: Connectome Workbench (https://www.humanconnectome.org/software/connectome-workbench)

## Instructions for Use

### Reproducibility

**Code:**

- Simulation studies: example case provided for N=500 and SNR = 1/0.2 (all other studies can also be reproduced by generating additional datasets with 'data_circle.R' or 'data_generation.R')
- UK Biobank application: Code to run ABCD study application where input data is concatenated into a NxM matrix from cifti files, where N represents sample size and M is the number of vertices, and save the output as a csv files and then in cifti files.


**Creation of figures and tables in the main paper and supplementary material:**

Paper:
- Fig. 1: gp_stgp_htgp_rtgp.R
- Fig. 2: data_generation.R
- Fig. 3: simstudy_surface.R
- Fig. 4: simstudy_surface.R
- Fig. 5: Created via Powerpoint.
- Fig. 6: (a) - (b) evaluation_RTGP.R, (c) - (d) GPR_Normal_GPR_Horseshoe.R, (e) - (f) BR_Normal_BR_Horseshoe.R, (g) - (h) Ridge_LASSO.R
- Fig. 7: (a) - (b) evaluation_RTGP.R, (c) - (d) GPR_Normal_GPR_Horseshoe.R
- Tab. 1,2,3: simstudy_surface.R
- Tab. 5: evaluation_RTGP.R

### Replication (Optional)

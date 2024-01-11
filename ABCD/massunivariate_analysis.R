# Script: Perform mass-univariate image-on-scalar analysis on ABCD data (X: intelligence scores, age, sex, ..., Y: test statistics of first-level
# fMRI analysis of emotional n-back 2 vs 0 back contrast) -> File is also used to derive subsampled true beta map for surface-based simulation study.

#############
### Paths ###
#############

# Path to Human Connectome workbench software.
wb_path = ''
# Path to folder where to store results in.
path_results = ''

# Create folder to store mass-univariate results in and change results path.
path_results = paste0(path_results, 'massunivariate/')
dir.create(file.path(path_results), showWarnings = FALSE)

# Create folder to store mass-univariate results with confounding variables in and change results path.
path_results = paste0(path_results, 'massunivariate/with_confounds/')
dir.create(file.path(path_results), showWarnings = FALSE)

# Path where ABCD data is stored in.
path_data = ''

# Path where an example of a cifti file is stored.
path_ABCD = ''

#################
### Libraries ###
#################

library(data.table)
library(MASS)
library(BayesGPfit)
library(parallel)
library(fields)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(truncnorm)
library(fastBayesReg)
library(STGP)
library(ggplot2)
library(ggpointdensity)
library(xtable)
library(bayestestR)

options(bitmapType='cairo-png')
library(ciftiTools)
ciftiTools.setOption('wb_path', wb_path)


#######################
### Left Hemisphere ###
#######################

# Data matrix with confounding variables.
X_confounds_train = fread(paste0(path_data, 'X_confounds_train.csv'))
X_confounds_train = data.matrix(X_confounds_train)[,2:ncol(X_confounds_train)]
# Data matrix with matrix(N x M, number of subjects x number of vertices) test statistics maps of emotional n-back task with 2 vs 0 back contrast
# (left hemisphere).
X_cortex_left_train = fread(paste0(path_data, 'X_left_cortex_train.csv'))
X_cortex_left_train = data.matrix(X_cortex_left_train)[,2:ncol(X_cortex_left_train)]
# Number of vertices.
M = ncol(X_cortex_left_train)
# Number of confounding variables.
P = ncol(X_confounds_train)
# Intelligence scores of training dataset.
y_train = as.numeric(unlist(fread(paste0(path_data, 'y_train.csv'))[,2]))

# Create matrices to store results in.
beta = matrix(NA, (P + 2), M)
beta_sd = matrix(NA, (P + 2), M)
p_value = matrix(NA, (P + 2), M)
test_statistic = matrix(NA, (P + 2), M)

# Perform a regression at each vertex j.
for(j in 1:M){
    print(j)
    model = lm(X_cortex_left_train[,j] ~ 1 + y_train + X_confounds_train)
    model = summary(model)
    beta[,j] = model$coefficients[,1]
    beta_sd[,j] = model$coefficients[,2]
    test_statistic[,j] = model$coefficients[,3]
    p_value[,j] = model$coefficients[,4]
}

# Determine binary significance based on p-values.
binary_significance = ifelse(p_value[2,] < 0.05, 1, 0)
# Determine binary significance based on FDR corrected p-values.
p_fdr = p.adjust(p_value, method = 'fdr', n = length(p_value))
p_fdr = matrix(p_fdr, (P+2), M)
binary_significance_corrected_fdr = ifelse(p_fdr[2,] < 0.05, 1, 0)
# Determine binary significance based on Bonferroni corrected p-values.
p_bonferroni = p.adjust(p_value, method = 'bonferroni', n = length(p_value))
p_bonferroni = matrix(p_bonferroni, (P+2), M)
binary_significance_corrected_bonferroni = ifelse(p_bonferroni[2,] < 0.05, 1, 0)

# Save results in list.
result_left_cortex = list()
result_left_cortex$beta = beta
result_left_cortex$beta_sd = beta_sd
result_left_cortex$test_statistic = test_statistic
result_left_cortex$p_value = p_value
result_left_cortex$p_fdr = p_fdr
result_left_cortex$p_bonferroni = p_bonferroni
result_left_cortex$binary_significance = binary_significance
result_left_cortex$binary_significance_corrected_fdr = binary_significance_corrected_fdr
result_left_cortex$binary_significance_corrected_bonferroni = binary_significance_corrected_bonferroni

save(result_left_cortex, file = paste0(path_results, "result_left_cortex.RData"))


########################
### Right Hemisphere ###
########################

# Data matrix with matrix(N x M, number of subjects x number of vertices) test statistics maps of emotional n-back task with 2 vs 0 back contrast
# (right hemisphere).
X_cortex_right_train = fread(paste0(path_data, 'X_right_cortex_train.csv'))
X_cortex_right_train = data.matrix(X_cortex_right_train)[,2:ncol(X_cortex_right_train)]
# Number of vertices.
M = ncol(X_cortex_right_train)

# Create matrices to store results in.
beta = matrix(NA, (P + 2), M)
beta_sd = matrix(NA, (P + 2), M)
p_value = matrix(NA, (P + 2), M)
test_statistic = matrix(NA, (P + 2), M)

# Perform a regression at each vertex j.
for(j in 1:M){
  print(j)
  model = lm(X_cortex_right_train[,j] ~ 1 + y_train + X_confounds_train)
  model = summary(model)
  beta[,j] = model$coefficients[,1]
  beta_sd[,j] = model$coefficients[,2]
  test_statistic[,j] = model$coefficients[,3]
  p_value[,j] = model$coefficients[,4]
}

# Determine binary significance based on p-values.
binary_significance = ifelse(p_value[2,] < 0.05, 1, 0)
# Determine binary significance based on FDR corrected p-values.
p_fdr = p.adjust(p_value, method = 'fdr', n = length(p_value))
p_fdr = matrix(p_fdr, (P+2), M)
binary_significance_corrected_fdr = ifelse(p_fdr[2,] < 0.05, 1, 0)
# Determine binary significance based on Bonferroni corrected p-values.
p_bonferroni = p.adjust(p_value, method = 'bonferroni', n = length(p_value))
p_bonferroni = matrix(p_bonferroni, (P+2), M)
binary_significance_corrected_bonferroni = ifelse(p_bonferroni[2,] < 0.05, 1, 0)

# Save results in list.
result_right_cortex = list()
result_right_cortex$beta = beta
result_right_cortex$beta_sd = beta_sd
result_right_cortex$test_statistic = test_statistic
result_right_cortex$p_value = p_value
result_right_cortex$p_fdr = p_fdr
result_right_cortex$p_bonferroni = p_bonferroni
result_right_cortex$binary_significance = binary_significance
result_right_cortex$binary_significance_corrected_fdr = binary_significance_corrected_fdr
result_right_cortex$binary_significance_corrected_bonferroni = binary_significance_corrected_bonferroni

save(result_right_cortex, file = paste0(path_results, "result_right_cortex.RData"))


# Read in results (left and right cortex)
load(paste0(path_results, 'result_left_cortex.RData'))
load(paste0(path_results, 'result_right_cortex.RData'))

# Save beta map of both hemispheres.
xii = read_xifti(paste0(path_ABCD, 'cope1.dtseries.nii'))
xii$data$cortex_left = matrix(result_left_cortex$beta[2,], length(result_left_cortex$beta[2,]), 1)
xii$data$cortex_right = matrix(result_right_cortex$beta[2,], length(result_right_cortex$beta[2,]), 1)

write_xifti(
  xii, 
  file.path(path_results, "beta_NM.dtseries.nii"), 
  file.path(path_results, "beta_NM_cortex_left.surf.gii"), file.path(path_results, "beta_NM_cortex_right.surf.gii")
)

# Save test statistics map of both hemispheres.
xii = read_xifti(paste0(path_ABCD, 'cope1.dtseries.nii'))
xii$data$cortex_left = matrix(result_left_cortex$test_statistic[2,], length(result_left_cortex$test_statistic[2,]), 1)
xii$data$cortex_right = matrix(result_right_cortex$test_statistic[2,], length(result_right_cortex$test_statistic[2,]), 1)

write_xifti(
  xii, 
  file.path(path_results, "test_statistic_NM.dtseries.nii"), 
  file.path(path_results, "test_statistic_NM_cortex_left.surf.gii"), file.path(path_results, "test_statistic_NM_cortex_right.surf.gii")
)

# Save binary significance map of both hemispheres.
xii = read_xifti(paste0(path_ABCD, 'cope1.dtseries.nii'))
xii$data$cortex_left = matrix(result_left_cortex$binary_significance, length(result_left_cortex$binary_significance), 1)
xii$data$cortex_right = matrix(result_right_cortex$binary_significance, length(result_right_cortex$binary_significance), 1)

write_xifti(
  xii, 
  file.path(path_results, "binary_significance_NM.dtseries.nii"), 
  file.path(path_results, "binary_significance_NM_cortex_left.surf.gii"), file.path(path_results, "binary_significance_NM_cortex_right.surf.gii")
)

# Save binary significance (FDR corrected) map of both hemispheres.
xii = read_xifti(paste0(path_ABCD, 'cope1.dtseries.nii'))
xii$data$cortex_left = matrix(result_left_cortex$binary_significance_corrected_fdr, length(result_left_cortex$binary_significance_corrected_fdr), 1)
xii$data$cortex_right = matrix(result_right_cortex$binary_significance_corrected_fdr, length(result_right_cortex$binary_significance_corrected_fdr), 1)

write_xifti(
  xii, 
  file.path(path_results, "binary_significance_corrected_fdr_NM.dtseries.nii"), 
  file.path(path_results, "binary_significance_corrected_fdr_NM_cortex_left.surf.gii"), file.path(path_results, "binary_significance_corrected_fdr_NM_cortex_right.surf.gii")
)

# Save binary significance (Bonferroni corrected) map of both hemispheres.
xii = read_xifti(paste0(path_ABCD, 'cope1.dtseries.nii'))
xii$data$cortex_left = matrix(result_left_cortex$binary_significance_corrected_bonferroni, length(result_left_cortex$binary_significance_corrected_bonferroni), 1)
xii$data$cortex_right = matrix(result_right_cortex$binary_significance_corrected_bonferroni, length(result_right_cortex$binary_significance_corrected_bonferroni), 1)

write_xifti(
  xii, 
  file.path(path_results, "binary_significance_corrected_bonferroni_NM.dtseries.nii"), 
  file.path(path_results, "binary_significance_corrected_bonferroni_NM_cortex_left.surf.gii"), file.path(path_results, "binary_significance_corrected_bonferroni_NM_cortex_right.surf.gii")
)


# Create true beta map for surface-based simulation study.
resampled_xii_fname = file.path(path_results, "beta_NM.dtseries.nii")

xii_2k = resample_cifti(
  cifti_fnames["dtseries"], resampled_xii_fname,
  resamp_res = 2000,
  write_dir = path_results
)

# Divide by 10.
xii_2k$data$cortex_left = xii_2k$data$cortex_left / 10
xii_2k$data$cortex_right = numeric(length(xii_2k$data$cortex_right))

write_xifti(
  xii_2k, 
  file.path(path_results, "beta_2k.dtseries.nii"), 
  file.path(path_results, "beta_2k.dtseries.nii_cortex_left.surf.gii"), file.path(path_results, "beta_2k.dtseries.nii_cortex_right.surf.gii")
)

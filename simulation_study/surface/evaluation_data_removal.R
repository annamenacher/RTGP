# Script: Plots, surface maps, and tables for simulation study evaluating the 10-fold CV.

#############
### Paths ###
#############

# Path to Human Connectome workbench software.
path_wb = ""
# Path to RTGP-VI model results with 100% of the data (without CV split)
path_100_percent = ""
# Path to RTGP-VI model results for each CV partition.
path_90_percent_core = path = ""
# Path to true beta map.
path_beta = ''

#################
### Libraries ###
#################

options(bitmapType='cairo-png')
library(ciftiTools)
ciftiTools.setOption('wb_path', path_wb)
library(data.table)
library(xtable)
library(ggplot2)

##################
### Evaluation ###
##################

# Load results of 100% data analysis of RTGP model.
load(paste0(path_100_percent, 'result_rtgp_vi_left_cortex.RData'))

# Determine highest ELBO value where algorithm converged.
iter = which(result$ELBO == max(result$ELBO))
# Binary significance values of 100% data analysis of RTGP.
alpha_thresholded_100 = result$alpha_thresholded[iter, ]

# Concatenate 10-fold CV results.
alpha_thresholded_90 = matrix(NA, 10, 1799)
for(sim in 1:10){
  # Read in results of CV partition.
  load(paste0(path_90_percent_core, 'sim', toString(sim), '/RTGP_VI_CV/result_rtgp_vi_left_cortex.RData'))
  # Determine highest ELBO value where algorithm converged.
  iter = which(result$ELBO == max(result$ELBO))
  # Binary significance results of CV partition.
  alpha_thresholded_90[sim,] = result$alpha_thresholded[iter, ]
}

# Selection stability across 10-fold CV
alpha_compare = apply(alpha_thresholded_90, 2, mean)

# Plot selection stability across 10-fold CV runs to check activation across CV partitions across vertices.
data = data.frame(x = 1:length(alpha_compare), y = alpha_compare)
ggplot(data, aes(x = x, y = y)) +
  geom_point() +
  labs(title = "Selection Stability",
       x = "Vertex",
       y = "Percentage of Selection") +
  theme_minimal()

# Save overlayed binary significance results of 10-fold CV.
xii = read_xifti(paste0(path_beta, 'beta_2k.dtseries.nii'))
xii$data$cortex_left = matrix(alpha_compare, length(alpha_compare), 1)
xii$data$cortex_right = matrix(0, 1803, 1)

write_xifti(
  xii, 
  file.path(path_90_percent_core, "alpha_thresholded_90_overlay.dtseries.nii"), 
  file.path(path_90_percent_core, "alpha_thresholded_90_overlay_left.surf.gii"), file.path(path_90_percent_core, "alpha_thresholded_90_overlay_cortex_right.surf.gii")
)

# Create matrices to save inference, parameter estimation, and prediction results for 10-fold CV. 
inference_df_all = matrix(NA, 10, 4)
parameters_df_all = matrix(NA, 10, 2)
prediction_df_all = matrix(NA, 10, 4)

for(sim in 1:10){

  path_out = paste0(path, 'sim', toString(sim), '/RTGP_VI_CV/')
  load(paste0(path_out, 'result_matrix.RData'))
  inference_df_all[sim,] = as.vector(as.numeric(unlist(result_matrix$matrix_inference[13,])))
  parameters_df_all[sim,] = as.vector(as.numeric(unlist(result_matrix$matrix_parameters[5,])))
  prediction_df_all[sim,] = as.vector(as.numeric(unlist(result_matrix$matrix_prediction[4,])))
  
}

# Read in results of 100% data analysis of RTGP model.
load(paste0(path_100_percent, "result_matrix.RData"))

# Inference results comparing 100% data analysis and averaged 90% data analysis in CV.
inference_df = matrix(NA, 2, 4)
inference_df[1,] = as.numeric(unlist(result_matrix$matrix_inference[7,]))
inference_df[2,] = apply(inference_df_all, 2, mean) 

# Parameter estimation results comparing 100% data analysis and averaged 90% data analysis in CV.
parameters_df = matrix(NA, 2, 2)
parameters_df[1,] = as.numeric(unlist(result_matrix$matrix_parameters[3,]))
parameters_df[2,] = apply(parameters_df_all, 2, mean) 

# Prediction results comparing 100% data analysis and averaged 90% data analysis in CV.
prediction_df = matrix(NA, 2, 4)
prediction_df[1,] = as.numeric(unlist(result_matrix$matrix_prediction[3,]))
prediction_df[2,] = apply(prediction_df_all, 2, mean) 

# Create tables for inference, prediction, and parameter estimation results.
row_names = c('100%', '90% (averaged)')
col_names_inference = c('TPR', 'TDR', 'FPR', 'FDR')
col_names_parameters = c('Bias', 'MSE')
col_names_prediction = c('R2 (train)', 'MSE (train)', 'R2 (test)', 'MSE (test)')

colnames(inference_df) = col_names_inference
rownames(inference_df) = row_names
colnames(parameters_df) = col_names_parameters
rownames(parameters_df) = row_names
colnames(prediction_df) = col_names_prediction
rownames(prediction_df) = row_names

# Inference results
xtable(inference_df, digits = 4)
print(xtable(inference_df, digits = 4), file=paste0(path, "inference.txt"))
write.csv(inference_df, paste0(path, 'matrix_inference.csv'))

# Parameter estimation results
xtable(parameters_df, digits = 4)
print(xtable(parameters_df, digits = 4), file=paste0(path, "parameters.txt"))
write.csv(parameters_df, paste0(path, 'matrix_parameters.csv'))

# Prediction results
xtable(prediction_df, digits = 4)
print(xtable(inference_df, digits = 4), file=paste0(path, "prediction.txt"))
write.csv(prediction_df, paste0(path, 'matrix_prediction.csv'))

# Save concatenated results.
write.csv(inference_df, paste0(path, 'inference.csv'))
write.csv(parameters_df, paste0(path, 'parameters.csv'))
write.csv(prediction_df, paste0(path, 'prediction.csv'))




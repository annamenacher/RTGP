# Script: Create tables with mean and standard error across results of 10 simulated datasets.

#############
### Paths ###
#############

path = '/well/nichols/users/fxo071/IOI-GP/sim_study/12.18.23/results/simstudy1/N500_L100_phi0.01_nu1.4/dataset'

path_out = '/well/nichols/users/fxo071/IOI-GP/sim_study/12.18.23/results/simstudy1/all/'

#################
### Libraries ###
#################

##################
### Evaluation ###
##################

# Create matrices to save mean and standard error across results of 10 simualted datasets.
matrix_inference = matrix(0, 10, 20)
matrix_parameter = matrix(0, 10, 14)
matrix_prediction = matrix(0, 10, 28)

for(i in 1: 10){
  
  # Load results.
  load(paste0(path, toString(i), '/all/result_matrix.RData'))
  # Store results.
  matrix_inference[i,] = matrix_inference[i,] + as.numeric(unlist(result_matrix$matrix_inference))
  matrix_parameter[i,] = matrix_parameter[i,] + as.numeric(unlist(result_matrix$matrix_parameter))
  matrix_prediction[i,] = matrix_prediction[i,] + as.numeric(unlist(result_matrix$matrix_prediction))
  
}

# Calculate mean of inference results across datasets and multiply by 100 for ease of interpretation.
matrix_inference_mean = round(apply(matrix_inference, 2, mean) * 100,2)
matrix_inference_mean = matrix(matrix_inference_mean, 5, 4)
colnames(matrix_inference_mean) = colnames(result_matrix$matrix_inference)
rownames(matrix_inference_mean) = rownames(result_matrix$matrix_inference)

# Calculate standard error of inference results across datasets and multiply by 100 for ease of interpretation.
matrix_inference_sd = round(apply(matrix_inference, 2, sd) * 100, 2)
matrix_inference_sd = matrix(matrix_inference_sd, 5, 4)
colnames(matrix_inference_sd) = colnames(result_matrix$matrix_inference)
rownames(matrix_inference_sd) = rownames(result_matrix$matrix_inference)

# Calculate mean of parameter estimates results across datasets and multiply by 100 for ease of interpretation.
matrix_parameter_mean = round(apply(matrix_parameter, 2, mean) * 100, 2)
matrix_parameter_mean = matrix(matrix_parameter_mean, 7, 2)
colnames(matrix_parameter_mean) = colnames(result_matrix$matrix_parameters)
rownames(matrix_parameter_mean) = rownames(result_matrix$matrix_parameters)

# Calculate standard error of parameter estimates results across datasets and multiply by 100 for ease of interpretation.
matrix_parameter_sd = round(apply(matrix_parameter, 2, sd) * 100, 2)
matrix_parameter_sd = matrix(matrix_parameter_sd, 7, 2)
colnames(matrix_parameter_sd) = colnames(result_matrix$matrix_parameters)
rownames(matrix_parameter_sd) = rownames(result_matrix$matrix_parameters)

# Calculate mean of prediction results across datasets and multiply by 100 for ease of interpretation.
matrix_prediction_mean = round(apply(matrix_prediction, 2, mean) * 100, 2)
matrix_prediction_mean = matrix(matrix_prediction_mean, 7, 4)
colnames(matrix_prediction_mean) = colnames(result_matrix$matrix_prediction)
rownames(matrix_prediction_mean) = rownames(result_matrix$matrix_prediction)

# Calculate standard error of prediction results across datasets and multiply by 100 for ease of interpretation.
matrix_prediction_sd = round(apply(matrix_prediction, 2, sd) * 100, 2)
matrix_prediction_sd = matrix(matrix_prediction_sd, 7, 4)
colnames(matrix_prediction_sd) = colnames(result_matrix$matrix_prediction)
rownames(matrix_prediction_sd) = rownames(result_matrix$matrix_prediction)

# Combine inference, parameter estimates, and prediction results.
matrix_inference_all = matrix(NA, 5, 8)
matrix_parameter_all = matrix(NA, 7, 4)
matrix_prediction_all = matrix(NA, 7, 8)

matrix_inference_all[1:5,1] = matrix_inference_mean[1:5,1]
matrix_inference_all[1:5,2] = matrix_inference_sd[1:5,1]
matrix_inference_all[1:5,3] = matrix_inference_mean[1:5,2]
matrix_inference_all[1:5,4] = matrix_inference_sd[1:5,2]
matrix_inference_all[1:5,5] = matrix_inference_mean[1:5,3]
matrix_inference_all[1:5,6] = matrix_inference_sd[1:5,3]
matrix_inference_all[1:5,7] = matrix_inference_mean[1:5,4]
matrix_inference_all[1:5,8] = matrix_inference_sd[1:5,4]




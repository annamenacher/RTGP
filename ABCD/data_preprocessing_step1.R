# Script: Preprocessing of data (step 1): Read in data from cifti files of emotional n-back task and split in left and right hemisphere data. 

#############
### Paths ###
#############

# Path to folder that contains Human Connectome workbench software.
wb_path = ''
# Path to folder where subject folders containing task fMRI data is located.
path_data = ''

#################
### Libraries ###
#################

library(ciftiTools)
ciftiTools.setOption('wb_path', wb_path)

######################
### Pre-processing ###
######################

# Read in subject IDs. 
ID_subj = read.csv("/well/nichols/shared/ABCD/cs-fMRI/CIFTI/nback_baseline.csv")
ID_subj = ID_subj$Subject
# Number of subjects with data available.
N = length(ID_subj)

# Create matrix (for left and right hemisphere) for beta coefficients from first level analysis of fMRI of emotional n-back task with 2 vs 0 back contrast. 
beta_left = matrix(NA, N, 29696)
beta_right = matrix(NA, N, 29716)

# Create matrix (for left and right hemisphere) for standard deviation of beta coefficients from first level analysis of fMRI of emotional n-back task 
# with 2 vs 0 back contrast. 
var_beta_left = matrix(NA, N, 29696)
var_beta_right = matrix(NA, N, 29716)

# Read in data from cifti file and store in matrix. 
for(i in 1:N){
  print(i)
  
  # Create current data path for specific subject i.
  path = paste0(path_data, toString(ID_subj[i]), "/baselineYear1Arm1/nback/task-nback_hp200_s2_level2.feat/GrayordinatesStats/cope11.feat/cope1.dtseries.nii")
  
  # Read in cifti file.
  xii = read_xifti(path)
  
  # Store COPE beta coefficients of left and right hemisphere in respective matrix. 
  beta_left[i,] = xii$data$cortex_left[,1]
  beta_right[i,] = xii$data$cortex_right[,1]
  
  # Create current data path for specific subject i.
  path = paste0(path_data, toString(ID_subj[i]), "/baselineYear1Arm1/nback/task-nback_hp200_s2_level2.feat/GrayordinatesStats/cope11.feat/varcope1.dtseries.nii")
 
  # Read in cifti file.
  xii = read_xifti(path)
  
  # Store standard deviation of beta coefficients of left and right hemisphere in respective matrix. 
  var_beta_left[i,] = xii$data$cortex_left[,1]
  var_beta_right[i,] = xii$data$cortex_right[,1]
  
}
 
# Save matrices of beta coefficients, standard deviation of beta coefficients and test statistics of left hemisphere data.
write.csv(beta_left, "/well/nichols/users/fxo071/IOI-GP/ABCD/data/test_statistics/beta_left_COPE.csv")
write.csv(var_beta_left, "/well/nichols/users/fxo071/IOI-GP/ABCD/data/test_statistics/beta_left_var_COPE.csv")
test_statistic_left = beta_left / sqrt(var_beta_left)
write.csv(test_statistic_left, "/well/nichols/users/fxo071/IOI-GP/ABCD/data/test_statistics/test_statistic_left.csv")

# Save matrices of beta coefficients, standard deviation of beta coefficients and test statistics of right hemisphere data.
write.csv(beta_right, "/well/nichols/users/fxo071/IOI-GP/ABCD/data/test_statistics/beta_right_COPE.csv")
write.csv(var_beta_right, "/well/nichols/users/fxo071/IOI-GP/ABCD/data/test_statistics/beta_right_var_COPE.csv")
test_statistic_right = beta_right / sqrt(var_beta_right)
write.csv(test_statistic_right, "/well/nichols/users/fxo071/IOI-GP/ABCD/data/test_statistics/test_statistic_right.csv")



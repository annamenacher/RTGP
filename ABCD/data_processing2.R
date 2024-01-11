# Script: Preprocessing of data (step 2): Create training, validation, and test dataset split of ABCD data.

#################
### Libraries ###
#################

library(data.table)
options(bitmapType='cairo-png')

#############
### Paths ###
#############

# Path to folder with test statistics matrices (derived from preprocessing step 1)
path_test_statistics = ''
# Path to folder with raw data (containing subject specific confounding variables, such as age, sex, ...)
path_raw_data = ''
# Path to file containing test statistic matrix of left hemisphere.
path_left_cortex = paste0(path_test_statistics, 'test_statistic_left.csv')
# Path to file containing test statistic matrix of right hemisphere.
path_right_cortex =  paste0(path_test_statistics, 'test_statistic_right.csv')
# Path to file containing raw data matrix.
path_raw_data = paste0(path_raw_data, 'raw_data.csv')
# Create folder to store unstandardised data in.
path_output_unstandardized =  paste0(path_test_statistics, 'not_centered/')
# Create folder to store centered data in.
path_output_standardized =  paste0(path_test_statistics, 'centered/')
# Create folder to store minimially processed data in (only remove NAs).
path_output_all =  paste0(path_test_statistics, 'all/') 

############
### Data ###
############

# Read in left hemisphere image data.
left_cortex = fread(path_left_cortex)
left_cortex = data.matrix(left_cortex)[,2:dim(left_cortex)[2]]

# Read in right hemisphere image data.
right_cortex = fread(path_right_cortex)
right_cortex = data.matrix(right_cortex)[,2:dim(right_cortex)[2]]

# Read in raw data with confounding variables and intelligence scores data.
raw_data = fread(path_raw_data)
raw_data = data.matrix(raw_data)[,2:dim(raw_data)[2]]

# Remove any subjects if there is any missing data in the imaging or raw data.
idx_remove = unique(c(which(is.na(raw_data), arr.ind = T)[,1], which(is.na(left_cortex), arr.ind = T)[,1], which(is.na(right_cortex), arr.ind = T)[,1]))
left_cortex = left_cortex[-idx_remove,]
right_cortex = right_cortex[-idx_remove,]
raw_data = raw_data[-idx_remove,]

# Pre-processing of raw data.
raw_data = raw_data[,2:23]

# Site scanner ID: categorical variable for the 21 scanner sites -> Scanner site 21: baseline scanner site.
site_id = matrix(NA, nrow(raw_data), 20)
site_names = rep('site.id.', 20)
for(i in 1:20){
  site_id[,i] = ifelse(raw_data[,4] == i, 1, 0)
  site_names[i] = paste0(site_names[i], toString(i))
}

# Remove not needed confounding variables.
raw_data = raw_data[,-4]
raw_data = raw_data[,-8]
raw_data = raw_data[,-8]
raw_data = raw_data[,-11]
raw_data = raw_data[,-c(14:18)]
raw_data_names = colnames(raw_data)
raw_data = cbind(raw_data, site_id)
colnames(raw_data) = c(raw_data_names, site_names)

# Output variable: intelligence scores
y = raw_data[,3]
# Matrix with confounding variables (parental income, parental marital status, age, sex, ...)
X_confounds = raw_data[,-3]

#####################
### Training data ###
#####################

# Set seed for reproducibility of training, validation, test split.
set.seed(1234)

# Number of subjects.
n = nrow(left_cortex)
# Number of vertices in left hemisphere.
M_left_cortex = ncol(left_cortex)
# Number of vertices in right hemisphere.
M_right_cortex = ncol(right_cortex)
# Number of confounding variables.
P = ncol(X_confounds)
  
# Number of subjects (training data)
n_train = round(n * 0.8, 0)
# Number of subjects (validation data)
n_val = round((n - n_train)/2, 0)
# Number of subjects (test data)
n_test = n - n_train - n_val
# Indices of training, validation, and test data split.
train_idx = sample(1:n, n_train)
val_idx = c(1:n)[-train_idx]
val_idx = val_idx[sample(1:(n_val+n_test), n_val)]
test_idx = c(1:n)[-c(train_idx, val_idx)]

# Save minimially processed data.
write.csv(left_cortex, paste0(path_output_all, 'left_cortex.csv'))
write.csv(right_cortex, paste0(path_output_all, 'right_cortex.csv'))
write.csv(raw_data, paste0(path_output_all, 'raw_data.csv'))

# Save training, validation, and test subject indices.
write.csv(train_idx, paste0(path_output_all, 'train_idx.csv'))
write.csv(test_idx, paste0(path_output_all, 'test_idx.csv'))
write.csv(val_idx, paste0(path_output_all, 'val_idx.csv'))

# Save unstandardised training dataset.
y_train = y[train_idx] 
X_confounds_train = X_confounds[train_idx,]
X_left_cortex_train = left_cortex[train_idx,]
X_right_cortex_train = right_cortex[train_idx,]

write.csv(y_train, paste0(path_output_unstandardized, 'y_train.csv'))
write.csv(X_confounds_train, paste0(path_output_unstandardized, 'X_confounds_train.csv'))
write.csv(X_left_cortex_train, paste0(path_output_unstandardized, 'X_left_cortex_train.csv'))
write.csv(X_right_cortex_train, paste0(path_output_unstandardized, 'X_right_cortex_train.csv'))

# Save mean centered training dataset.
y_train_standarized = (y_train - mean(y_train)) 

interview_age = (X_confounds[train_idx,2] - mean(X_confounds[train_idx,2]))
X_confounds_train_standardized = X_confounds_train
X_confounds_train_standardized[,2] = interview_age

mean_left_cortex_train = as.vector(apply(X_left_cortex_train, 2, mean))
X_left_cortex_train_standardized = sweep(X_left_cortex_train, 2, mean_left_cortex_train, "-")

mean_right_cortex_train = as.vector(apply(X_right_cortex_train, 2, mean))
X_right_cortex_train_standardized = sweep(X_right_cortex_train, 2, mean_right_cortex_train, "-")

write.csv(y_train_standarized, paste0(path_output_standardized, 'y_train_centered.csv'))
write.csv(mean(y_train), paste0(path_output_standardized, 'y_train_mean.csv'))

write.csv(X_confounds_train_standardized, paste0(path_output_standardized, 'X_confounds_train_centered.csv'))
write.csv(mean(X_confounds_train[,2]), paste0(path_output_standardized, 'interview_age_train_mean.csv'))

write.csv(X_left_cortex_train_standardized, paste0(path_output_standardized, 'X_left_cortex_train_centered.csv'))
write.csv(mean_left_cortex_train, paste0(path_output_standardized, 'mean_left_cortex_train.csv'))

write.csv(X_right_cortex_train_standardized, paste0(path_output_standardized, 'X_right_cortex_train_centered.csv'))
write.csv(mean_right_cortex_train, paste0(path_output_standardized, 'mean_right_cortex_train.csv'))

# Save unstandardised test dataset.
y_test = y[test_idx] 
X_confounds_test = X_confounds[test_idx,]
X_left_cortex_test = left_cortex[test_idx,]
X_right_cortex_test = right_cortex[test_idx,]

write.csv(y_test, paste0(path_output_unstandardized, 'y_test.csv'))
write.csv(X_confounds_test, paste0(path_output_unstandardized, 'X_confounds_test.csv'))
write.csv(X_left_cortex_test, paste0(path_output_unstandardized, 'X_left_cortex_test.csv'))
write.csv(X_right_cortex_test, paste0(path_output_unstandardized, 'X_right_cortex_test.csv'))

# Save mean centered test dataset.
y_test_standarized = (y_test - mean(y_train)) 

interview_age = (X_confounds[test_idx,2] - mean(X_confounds[train_idx,2])) 
X_confounds_test_standardized = X_confounds_test
X_confounds_test_standardized[,2] = interview_age

X_left_cortex_test_standardized = sweep(X_left_cortex_test, 2, mean_left_cortex_train, "-")

X_right_cortex_test_standardized = sweep(X_right_cortex_test, 2, mean_right_cortex_train, "-")

write.csv(y_test_standarized, paste0(path_output_standardized, 'y_test_centered.csv'))
write.csv(X_confounds_test_standardized, paste0(path_output_standardized, 'X_confounds_test_centered.csv'))
write.csv(X_left_cortex_test_standardized, paste0(path_output_standardized, 'X_left_cortex_test_centered.csv'))
write.csv(X_right_cortex_test_standardized, paste0(path_output_standardized, 'X_right_cortex_test_centered.csv'))

# Save unstandardised validation dataset.
y_val = y[val_idx] 
X_confounds_val = X_confounds[val_idx,]
X_left_cortex_val = left_cortex[val_idx,]
X_right_cortex_val = right_cortex[val_idx,]

write.csv(y_val, paste0(path_output_unstandardized, 'y_val.csv'))
write.csv(X_confounds_val, paste0(path_output_unstandardized, 'X_confounds_val.csv'))
write.csv(X_left_cortex_val, paste0(path_output_unstandardized, 'X_left_cortex_val.csv'))
write.csv(X_right_cortex_val, paste0(path_output_unstandardized, 'X_right_cortex_val.csv'))

# Save mean centered validation dataset.
y_val_standarized = (y_val - mean(y_train)) 

interview_age = (X_confounds[val_idx,2] - mean(X_confounds[train_idx,2])) 
X_confounds_val_standardized = X_confounds_val
X_confounds_val_standardized[,2] = interview_age

X_left_cortex_val_standardized = sweep(X_left_cortex_val, 2, mean_left_cortex_train, "-")

X_right_cortex_val_standardized = sweep(X_right_cortex_val, 2, mean_right_cortex_train, "-")

write.csv(y_val_standarized, paste0(path_output_standardized, 'y_val_centered.csv'))
write.csv(X_confounds_val_standardized, paste0(path_output_standardized, 'X_confounds_val_centered.csv'))
write.csv(X_left_cortex_val_standardized, paste0(path_output_standardized, 'X_left_cortex_val_centered.csv'))
write.csv(X_right_cortex_val_standardized, paste0(path_output_standardized, 'X_right_cortex_val_centered.csv'))



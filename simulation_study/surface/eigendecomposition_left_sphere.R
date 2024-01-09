# Script: Perform eigendecomposition of kernel matrix for left hemisphere of the brain (for various kernel hyperparameters)

sim = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

#################
### Libraries ###
#################

library(RSpectra)
library(data.table)

#################
### Constants ###
#################

# Kernel hyperparameter: Bandwidth
phi_list = c(0.01, 0.05, seq(0.1, 1, length.out = 10))
phi = phi_list[sim]
# Kernel hyperparameter: Exponent
nu_list = c(1.2, 1.4, 1.6, 1.8, 2)
nu = nu_list[sim]
# Number of basis functions
L = 100

#############
### Paths ###
#############

# Path to folder where to save the results of the eigendecomposition and where the distance matrix is stored.
path = ''
path_out = paste0(path, 'eigendecomposition_L', toString(L), '_phi', toString(phi), '_nu', toString(nu), '/')
dir.create(file.path(path_out), showWarnings = FALSE)

# Distance matrix 
distance_matrix = data.matrix(fread(paste0(path, 'LeftDenseGeodesicDistmat.txt')))

if(dim(distance_matrix)[1] != dim(distance_matrix)[2]){
  distance_matrix = distance_matrix[,2:dim(distance_matrix)[2]]
}

# Delete locations in the distance matrix that are not actually in the data (compared to the template)
mask = as.numeric(unlist(fread('/well/nichols/users/fxo071/IOI-GP/sim_study/11.21.23/data/idx_left.csv')[,2]))
distance_matrix = distance_matrix[mask, mask]

# Two parameter exponential radial basis kernel function
covariance_kernel = function(x, s, phi, nu){
  y = s * exp(-phi*x^nu)
  return(y)
}

# Calculate kernel matrix
distance_matrix = covariance_kernel(distance_matrix, s = 1, phi = phi, nu = nu)

# Eigendecomposition of kernel matrix
eigendecomp = eigs_sym(A = distance_matrix, k = L, which = "LM") 

rm(distance_matrix)

# Save results of eigendecomposition
write.csv(eigendecomp$values, paste0(path_out, 'eigenvalues_left_cortex.csv'))
write.csv(eigendecomp$vectors, paste0(path_out, 'eigenvectors_left_cortex.csv'))
save(eigendecomp, file = paste0(path_out, 'eigendecomposition_left_cortex.RData'))

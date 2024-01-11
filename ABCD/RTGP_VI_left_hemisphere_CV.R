# Script: Perform model estimation of RTGP with variational inference (for left hemisphere, 10-fold cross validation)

sim = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Kernel hyperparameter settings
# Number of basis functions
L = 800
# Bandwidth parameter
phi = 0.2
# Exponent parameter
nu = 1.2

# Remove 10% of the training data (10-fold cross validation)
percent_removal = 10
# Approximately 300 subjects in the test dataset.
n_remove = 300

#############
### Paths ###
#############

# Path to folder containing Human connectome workbench software.
wb_path = ''

# Path to folder containing the data.
path_data = ""

# Path to folder containing the eigendecompositions
path_eigendecomp = ""

# Path to folder where to store results in.
path_out = path = ""
dir.create(file.path(path), showWarnings = FALSE)

path_out = path = paste0(path_out, toString(sim), "/")
dir.create(file.path(path), showWarnings = FALSE)

# Path to folder where initial parameters are stored in (use GPR + Normal estimates) 
path_init = ""

#################
### Libraries ###
#################

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
library(data.table)
library(ciftiTools)
options(bitmapType='cairo-png')
ciftiTools.setOption('wb_path', wb_path)

############
### Data ###
############

# Read in data matrix with confounding variables (training data).
X_confounds = fread(paste0(path_data, 'X_confounds_train.csv'))
X_confounds = data.matrix(X_confounds)[,2:ncol(X_confounds)]
# Read in imaging data of left hemisphere (training data).
X = fread(paste0(path_data, 'X_left_cortex_train.csv'))
X = data.matrix(X)[,2:ncol(X)]
# Read in output vector of intelligence scores (training data)
y = as.numeric(unlist(fread(paste0(path_data, 'y_train.csv'))[,2]))

# Perform cross validation split into 10 partitions and identify indices.
idx = c(((sim - 1)*n_remove + 1) : (n_remove*sim))

# Split training data into 10-fold CV partition.
X = X[-idx, ]
X_confounds = X_confounds[-idx, ]
y = y[-idx]

# Save the index of test data of each partition.
write.csv(idx, paste0(path_out, 'idx.csv'))

# Read in eigenfunctions and eigenvalues of eigendecomposition.
Psi = fread(paste0(path_eigendecomp, '/eigenvectors_left_cortex.csv'))
Psi = t(data.matrix(Psi)[,2:ncol(Psi)])
lambda_l = as.numeric(unlist(fread(paste0(path_eigendecomp, 'eigenvalues_left_cortex.csv'))[,2]))

# Number of vertices
M = p = ncol(X)
# Number of confounding variables
P = ncol(X_confounds)
# Number of subjects in training dataset
n_train = nrow(X_confounds)

# Read in list with results of GPR + Normal model and use the estimated parameters as initialisation.
load(paste0(path_init, "/result_left_cortex.RData"))
beta0_NM = result_left_cortex$nm_fit_SOI$post_mean$betacoef[1,]
beta_confound_NM = result_left_cortex$nm_fit_SOI$post_mean$betacoef[2:(P+1),]
beta_cortex_NM = result_left_cortex$beta_fit$NM
theta_NM = result_left_cortex$theta_fit$NM
sigma_eps2_NM = result_left_cortex$nm_fit_SOI$post_mean$sigma2_eps


###############
### RTGP-VI ###
###############

# Enforce maximum number of iterations in case algorithm does not converge.
n_iter = 1000
# Convergence criterion: If difference of ELBO values is lower than this threshold stop algorithm.
eps = 0.00001

# Standard normal density function.
phi_function = function(x){
  x = (1/sqrt(2*pi)) * exp((-0.5)*x^2)
  return(x)
}

# Log-likelihood function of RTGP model.
log_likelihood = function(y, Psi, lambda_l, theta, beta0, X, sigma_epsilon, alpha, threshold){
  beta = t(Psi) %*% diag(as.vector(sqrt(lambda_l)), L) %*% theta
  beta_thresh = matrix(beta * ifelse(abs(alpha) > threshold, 1, 0), p, 1)
  mu = beta0 + X%*%beta_thresh + X_confounds%*%beta_confound
  ll = sum(dnorm(y, mu, sigma_epsilon, log=TRUE))
  return(ll)
}

# Joint density function of RTGP model.
joint_density = function(){
  beta = t(Psi) %*% diag(as.vector(sqrt(lambda_l)), L) %*% theta
  beta_thresh = matrix(beta * ifelse(abs(alpha) > threshold, 1, 0), p, 1)
  mu = beta0 + X%*%beta_thresh + X_confounds %*% beta_confound
  joint = sum(dnorm(y, mu, sigma_epsilon, log=TRUE)) + dnorm(beta0, 0, sigma_beta0, log = T) + 
    sum(dnorm(theta, 0, sigma_beta, log = T)) + 
    sum(dnorm(alpha, beta, sigma_alpha, log = T)) +
    (-1)*(a_beta + 1)*log(sigma_beta^2) - ((b_beta) / sigma_beta^2)  +
    (-1)*(a_alpha + 1)*log(sigma_alpha^2) - ((b_alpha) / sigma_alpha^2)  +
    (-1)*(a_epsilon + 1)*log(sigma_epsilon^2) - ((b_epsilon) / sigma_epsilon^2) 
  return(joint)
}

# Hyperparameters
# Hyperparameters of prior on variances sigma_alpha, sigma_beta, sigma_epsilon.
a_beta = b_beta = 0.001
a_alpha = b_alpha = 0.001
a_epsilon = b_epsilon = 0.001
# Lower and upper bound of discrete Uniform prior on threshold lambda.
t_min = max(beta_cortex_NM) * 0.25
t_max =  max(beta_cortex_NM) * 0.9
# Threshold options of discrete Uniform prior on threshold lambda.
threshold_options = seq(t_min , t_max, length.out = 10)
# Hyperparamter of prior on intercept beta0
sigma_beta0 = 10
# Hyperparameter of prior on confounding variables
sigma_beta_confound = 10

# Initialisation of parameters
theta = matrix(theta_NM,L,1)
beta0 =  beta0_NM
beta = beta_cortex_NM
beta_confound = beta_confound_NM
alpha = beta_cortex_NM
threshold = (t_min + t_max) / 2
sigma_epsilon = sigma_epsilon_prior = sqrt(sigma_eps2_NM)
sigma_alpha = sigma_alpha_prior = 0.1
sigma_beta = sigma_beta_prior = 1

# Initialise other quantities needed for VI updates.
# Inverse of variance parameters.
exp_1_sigma_epsilon2 = 1/sigma_epsilon^2
exp_1_sigma_alpha2 = 1/sigma_alpha^2
exp_1_sigma_beta2 = 1/sigma_beta^2
# Posterior variance of variational density on intercept.
Sigma_beta0 = 1/(n_train*exp_1_sigma_epsilon2 + (1/sigma_beta0^2))
# Posterior variance of variational density on theta.
Sigma_theta = solve(exp_1_sigma_epsilon2 * (t(X%*%diag(ifelse(abs(as.vector(alpha))>threshold,1,0),p)%*%t(Psi)%*%diag(sqrt(lambda_l),L))%*%X%*%diag(ifelse(abs(as.vector(alpha))>threshold,1,0),p)%*%t(Psi)%*%diag(sqrt(lambda_l),L) + diag(exp_1_sigma_beta2,L) + exp_1_sigma_alpha2*t(t(Psi)%*%diag(sqrt(lambda_l),L))%*%t(Psi)%*%diag(sqrt(lambda_l),L)))
# Expected value of indictor function of absolute value of alpha greater than threshold lambda.
exp_I_abs_alpha_greater_lamda = ifelse(abs(alpha)>threshold, 1, 0)
# Expected value of beta^2
exp_beta2 = diag((t(Psi) %*% diag(sqrt(lambda_l),L)) %*% (Sigma_theta + theta %*% t(theta)) %*% t(t(Psi) %*% diag(sqrt(lambda_l),L)))
# Reweighting of weights on mixture of truncated normal distributions on alpha.
c = exp((-0.5) * exp_1_sigma_epsilon2 * 2 * t(t(y - beta0 - X_confounds%*%beta_confound) %*% (-X)) * beta  - 0.5 * exp_1_sigma_epsilon2 * colSums(X^2) * exp_beta2-0.5 * exp_1_sigma_epsilon2 * colSums(sweep(X, 2, beta, "*")*(sweep(-sweep(X, 2, (beta * exp_I_abs_alpha_greater_lamda), "*") , 1, X %*% (beta * exp_I_abs_alpha_greater_lamda), "+"))))
# Weights of mixture distribution on alpha.
w = matrix(NA, 3, p)
w[1,] = as.vector(pnorm(-(threshold + beta)*sqrt(exp_1_sigma_alpha2))) * c
w[2,] = as.vector(pnorm((threshold - beta)*sqrt(exp_1_sigma_alpha2)) - pnorm(-(threshold + beta)*sqrt(exp_1_sigma_alpha2)))
w[3,] = as.vector(1 - pnorm((threshold - beta)*sqrt(exp_1_sigma_alpha2))) * c
w_sum = apply(w, 2, sum)
w = sweep(w, 2, w_sum, FUN="/")
# Expected value of alpha.
exp_alpha = matrix(NA, 3, p)
exp_alpha[1,] = beta - sigma_alpha * (phi_function((-threshold - beta) / sigma_alpha) / pnorm((-threshold - beta) / sigma_alpha))
exp_alpha[2,] = beta - sigma_alpha * (phi_function((threshold - beta) / sigma_alpha) - phi_function((-threshold - beta) / sigma_alpha)) / (pnorm((threshold - beta) / sigma_alpha) - pnorm((-threshold - beta) / sigma_alpha))
exp_alpha[3,] = beta + sigma_alpha * phi_function((threshold - beta) / sigma_alpha) / (1 - pnorm((threshold - beta) / sigma_alpha))

# Initialise difference of ELBOs.
diff = 100
# Initialise iteration counter of algorithm.
iter = 0
# Initialise ELBO value.
ELBO = -100000000000000000

# Create matrices to store results in across iterations. 
log_likelihood_chain = numeric(n_iter)
joint_chain = numeric(n_iter)
diff_chain = numeric(n_iter)
threshold_chain = numeric(n_iter)
ELBO_chain = numeric(n_iter)
beta0_chain = numeric(n_iter)
beta_confound_chain = matrix(NA, n_iter, P)
theta_chain = matrix(NA, n_iter, L)
alpha_chain = matrix(NA, n_iter, M)
beta_chain = matrix(NA, n_iter, M)
beta_threshold_chain = matrix(NA, n_iter, M)
alpha_threshold_chain = matrix(NA, n_iter, M)
sigma_epsilon_chain = numeric(n_iter)
sigma_alpha_chain = numeric(n_iter)
sigma_beta_chain = numeric(n_iter)

# Start time of analysis.
time.start = Sys.time()

while((diff > eps) && (iter != n_iter)){
  tryCatch({
    
    # Update iteration counter.
    iter = iter + 1
    print(iter)
    
    # Save results of previous iterations.
    beta_old = beta
    alpha_old = alpha
    ELBO_old = ELBO
    
    # Calculate log-likelihood and joint density value of RTGP model.
    ll = log_likelihood(y, Psi, lambda_l, theta, beta0, X, sigma_epsilon, alpha, threshold)
    joint = joint_density()
    
    # Update beta0 
    exp_I_abs_alpha_greater_lamda = w[1,] + w[3,]
    Sigma_beta0 = 1/(n_train*exp_1_sigma_epsilon2 + (1/(sigma_beta0^2)))
    beta0 = as.vector( Sigma_beta0 * (exp_1_sigma_epsilon2 * sum(y - X_confounds %*% beta_confound - X%*%(exp_I_abs_alpha_greater_lamda * beta))))
    
    # Update beta_confound 
    Sigma_beta_confound = solve(diag(1/sigma_beta_confound^2, P) + t(X_confounds) %*% X_confounds)
    beta_confound = Sigma_beta_confound %*% (t(X_confounds)%*%(y - beta0 - X%*%(exp_I_abs_alpha_greater_lamda * beta)))
      
    # Update theta
    Sigma_theta = solve(exp_1_sigma_epsilon2 * t(X %*% diag(as.vector(exp_I_abs_alpha_greater_lamda), p) %*% t(Psi) %*% diag(as.vector(sqrt(lambda_l)), L)) %*% (X %*% diag(as.vector(exp_I_abs_alpha_greater_lamda), p) %*% t(Psi) %*% diag(as.vector(sqrt(lambda_l)), L)) + diag(exp_1_sigma_beta2, L) + exp_1_sigma_alpha2*t(t(Psi)%*%diag(as.vector(sqrt(lambda_l)), L)) %*% (t(Psi)%*%diag(as.vector(sqrt(lambda_l)), L)))
    theta = Sigma_theta %*% t(exp_1_sigma_epsilon2 * t(y - beta0 - X_confounds%*%beta_confound) %*% (X%*%diag(as.vector(exp_I_abs_alpha_greater_lamda),p)%*%t(Psi)%*%diag(as.vector(sqrt(lambda_l)),L)) + exp_1_sigma_alpha2 * t(alpha)%*%t(Psi)%*%diag(as.vector(sqrt(lambda_l)),L))  
    beta = t(Psi)%*%(sqrt(lambda_l) * theta)

    # Update alpha
    exp_beta2 = diag((t(Psi) %*% diag(sqrt(lambda_l),L)) %*% (Sigma_theta + theta %*% t(theta)) %*% t(t(Psi) %*% diag(sqrt(lambda_l),L)))
    c = exp((-0.5) * exp_1_sigma_epsilon2 * 2 * t(t(y - beta0 - X_confounds%*%beta_confound) %*% (-X)) * beta  - 0.5 * exp_1_sigma_epsilon2 * colSums(X^2) * exp_beta2-0.5 * exp_1_sigma_epsilon2 * colSums(sweep(X, 2, beta, "*")*(sweep(-sweep(X, 2, (beta * exp_I_abs_alpha_greater_lamda), "*") , 1, X %*% (beta * exp_I_abs_alpha_greater_lamda), "+"))))

    w = matrix(NA, 3, p)
    w[1,] = as.vector(pnorm(-(threshold + beta)*sqrt(exp_1_sigma_alpha2))) * c
    w[2,] = as.vector(pnorm((threshold - beta)*sqrt(exp_1_sigma_alpha2)) - pnorm(-(threshold + beta)*sqrt(exp_1_sigma_alpha2)))
    w[3,] = as.vector(1 - pnorm((threshold - beta)*sqrt(exp_1_sigma_alpha2))) * c
    w_sum = apply(w, 2, sum)
    w = sweep(w, 2, w_sum, FUN="/")

    exp_alpha = matrix(NA, 3, p)
    exp_alpha[1,] = beta - sigma_alpha * (phi_function((-threshold - beta) / sigma_alpha) / (pnorm((-threshold - beta) / sigma_alpha) + 10^(-10)))
    exp_alpha[2,] = beta - sigma_alpha * (phi_function((threshold - beta) / sigma_alpha) - phi_function((-threshold - beta) / sigma_alpha)) / (pnorm((threshold - beta) / sigma_alpha) - pnorm((-threshold - beta) / sigma_alpha) + 10^(-10))
    exp_alpha[3,] = beta + sigma_alpha * phi_function((threshold - beta) / sigma_alpha) / (1 - pnorm((threshold - beta) / sigma_alpha) + 10^(-10))


    alpha = matrix(w[1,] * exp_alpha[1,] + w[2,] * exp_alpha[2,] + w[3,] * exp_alpha[3,], p, 1)    #alpha = exp_alpha


    # Update threshold
    threshold_t = numeric(length(threshold_options))
    for(k in 1:length(threshold_options)){
      threshold_t[k] = (-0.5) * exp_1_sigma_epsilon2 * sum((y - beta0 - X_confounds%*%beta_confound - X%*%(beta * ifelse(abs(alpha) > threshold_options[k], 1, 0)))^2)
    }

    logProb = threshold_t - max(threshold_t)
    Prob = exp(logProb)/sum(exp(logProb))
    threshold = sum(Prob * threshold_options)

    # Update sigma_epsilon
    w = matrix(NA, 3, p)
    w[1,] = as.vector(pnorm(-(threshold + beta)*sqrt(exp_1_sigma_alpha2))) * c
    w[2,] = as.vector(pnorm((threshold - beta)*sqrt(exp_1_sigma_alpha2)) - pnorm(-(threshold + beta)*sqrt(exp_1_sigma_alpha2)))
    w[3,] = as.vector(1 - pnorm((threshold - beta)*sqrt(exp_1_sigma_alpha2))) * c
    w_sum = apply(w, 2, sum)
    w = sweep(w, 2, w_sum, FUN="/")
    exp_I_abs_alpha_greater_lamda = w[1,] + w[3,]
    var_beta_alpha = exp_beta2 * exp_I_abs_alpha_greater_lamda - beta^2 * exp_I_abs_alpha_greater_lamda^2
    exp_beta_alpha_X = sum(sweep(X^2, 2, var_beta_alpha, "*")) + sum(diag((X%*%(beta * exp_I_abs_alpha_greater_lamda)) %*% t(X%*%(beta * exp_I_abs_alpha_greater_lamda))))
    
    # Update sigma_epsilon
    exp_1_sigma_epsilon2 = as.vector((n_train + a_epsilon) / (b_epsilon + 0.5 * sum((y - beta0 - X_confounds%*%beta_confound - X%*%(beta * exp_I_abs_alpha_greater_lamda))^2)))
    sigma_epsilon = as.vector(sqrt(1/exp_1_sigma_epsilon2))
    
    # Update sigma_beta
    exp_1_sigma_beta2 = as.vector((L/2 + a_beta) / (0.5 * sum(diag(Sigma_theta) + theta^2) + b_beta))
    sigma_beta = as.vector(sqrt(1/exp_1_sigma_beta2))
    
    # Update sigma_alpha
    var_alpha = matrix(NA, 3, p)
    var_alpha[1,] = sigma_alpha^2 * (1 - ((-threshold - beta)/sigma_alpha) * phi_function(((-threshold - beta)/sigma_alpha)) / (pnorm(((-threshold - beta)/sigma_alpha)) + 10^(-10)) - (phi_function(((-threshold - beta)/sigma_alpha)) / (pnorm(((-threshold - beta)/sigma_alpha)) + 10^(-10)))^2)
    var_alpha[2,] = sigma_alpha^2 * (1 - (((threshold - beta)/sigma_alpha) * phi_function(((threshold - beta)/sigma_alpha)) - ((-threshold - beta)/sigma_alpha) * phi_function(((-threshold - beta)/sigma_alpha))) / (pnorm(((threshold - beta)/sigma_alpha)) - pnorm(((-threshold - beta)/sigma_alpha)) + 10^(-10)) - ((phi_function(((threshold - beta)/sigma_alpha)) - phi_function(((-threshold - beta)/sigma_alpha))) / (pnorm(((threshold - beta)/sigma_alpha)) - pnorm(((-threshold - beta)/sigma_alpha)) + 10^(-10)))^2)
    var_alpha[3,] = sigma_alpha^2 * (1 + (((threshold - beta)/sigma_alpha)) * phi_function(((threshold - beta)/sigma_alpha)) / (1 - pnorm(((threshold - beta)/sigma_alpha)) + 10^(-10)) - (phi_function(((threshold - beta)/sigma_alpha)) / (1 - pnorm(((threshold - beta)/sigma_alpha)) + 10^(-10)))^2)

    exp_alpha2 = w[1,] * (var_alpha[1,] + exp_alpha[1,]^2) + w[2,] * (var_alpha[2,] + exp_alpha[2,]^2) + w[3,] * (var_alpha[3,] + exp_alpha[3,]^2)
    exp_1_sigma_alpha2 = (p/2 + a_alpha) / (b_alpha + 0.5 * sum((alpha - beta)^2))
    sigma_alpha = as.vector(sqrt(1/exp_1_sigma_alpha2))

    # Calculate thresholded alpha and beta values 
    beta_thresholded =  as.vector( beta * ifelse(exp_I_abs_alpha_greater_lamda > 0.5, 1, 0))
    alpha_thresholded = as.vector( ifelse(exp_I_abs_alpha_greater_lamda > 0.5, 1, 0) )

    # Save log-likelihood, joint and threshold values in vectors.
    log_likelihood_chain[(iter)] = ll
    joint_chain[iter] = joint
    threshold_chain[iter] = threshold
    
    # Calculate ELBO
    exp_ln_sigma_alpha2 = log(0.5 * sum((alpha-beta)^2) + b_alpha) - digamma(p/2 + a_alpha)
    exp_ln_sigma_epsilon2 = log(0.5 * (t(y)%*%y + n_train * (Sigma_beta0 + beta0^2) + exp_beta_alpha_X - 2 * t(y-beta0)%*%(X%*%(beta * exp_I_abs_alpha_greater_lamda)) - sum(2* beta0 * y)) + b_epsilon) - digamma(n_train/2 + a_epsilon)
    exp_ln_sigma_beta2 = log(0.5 * sum(diag(Sigma_theta) + theta^2) + b_beta) - digamma(L/2 + a_beta)
    part = log(pnorm((-threshold - beta)*sqrt(exp_1_sigma_alpha2)) + 10^(-10))
    part[part == -Inf] = 0
    part[part == Inf] = 0
    alpha_part1 = sum(w[1,] * part)
    part = -(-threshold - beta)*sqrt(exp_1_sigma_alpha2) * phi_function((-threshold - beta)*sqrt(exp_1_sigma_alpha2)) / (2*pnorm((-threshold - beta)*sqrt(exp_1_sigma_alpha2)) + 10^(-10))
    part[part == -Inf] = 0
    part[part == Inf] = 0
    alpha_part1 = alpha_part1 + sum(w[1,] * part)
    part = log(pnorm((threshold - beta)*sqrt(exp_1_sigma_alpha2)) - pnorm((-threshold-beta)*sqrt(exp_1_sigma_alpha2)) + 10^(-10))
    part[part == -Inf] = 0
    part[part == Inf] = 0
    alpha_part2 = sum(w[2,] * part)
    part = ((-threshold - beta)*sqrt(exp_1_sigma_alpha2) * phi_function((-threshold - beta)*sqrt(exp_1_sigma_alpha2)) - (threshold - beta)*sqrt(exp_1_sigma_alpha2)*phi_function((threshold - beta)*sqrt(exp_1_sigma_alpha2))) / (2 * (pnorm((threshold - beta)*sqrt(exp_1_sigma_alpha2)) - pnorm((-threshold - beta)*sqrt(exp_1_sigma_alpha2))) + 10^(-10))
    part[part == -Inf] = 0
    part[part == Inf] = 0
    alpha_part2 = alpha_part2 + sum(w[2,] * part)
    part = log(1 - pnorm((threshold-beta)*sqrt(exp_1_sigma_alpha2)) + 10^(-10))
    part[part == -Inf] = 0
    part[part == Inf] = 0
    alpha_part3 = sum(w[3,] * part)
    part = ((threshold-beta)*sqrt(exp_1_sigma_alpha2)*phi_function((threshold-beta)*sqrt(exp_1_sigma_alpha2))) / (2*(1-pnorm((threshold-beta)*sqrt(exp_1_sigma_alpha2))) + 10^(-10))
    part[part == -Inf] = 0
    part[part == Inf] = 0
    alpha_part3 = alpha_part3 + sum(w[3,] * part)

    ELBO = - (n_train/2) * log(2*pi) - (n_train/2) * exp_ln_sigma_epsilon2 - 0.5 * exp_1_sigma_epsilon2 * sum((y - beta0 - X_confounds%*%beta_confound - X%*%beta_thresholded )^2) -
      0.5 * log(2*pi) - 0.5 * log(sigma_beta0^2) - (1/sigma_beta0^2) * (Sigma_beta0 + beta0^2) - 0.5 * P * log(2*pi) - 0.5 * P * log(sigma_beta_confound^2) - (1/sigma_beta_confound^2) * sum(diag(Sigma_beta_confound + beta_confound%*%t(beta_confound))) -
      (L/2) * log(2*pi) - (L/2) * exp_ln_sigma_beta2 -
      0.5 * exp_1_sigma_beta2 * sum(diag(Sigma_theta) + theta^2) - (p/2) * log(2*pi) - (p/2) * exp_ln_sigma_alpha2 - 0.5 * exp_1_sigma_alpha2 * sum((alpha-beta)^2) +
      a_epsilon * log(b_epsilon) + (a_epsilon + 1) * exp_ln_sigma_epsilon2 - b_epsilon * exp_1_sigma_epsilon2 + a_beta * log(b_beta) +
      (a_beta + 1) * exp_ln_sigma_beta2 - b_beta * exp_1_sigma_beta2 + a_alpha * log(b_alpha) + (a_alpha + 1) * exp_ln_sigma_alpha2 -
      b_alpha * exp_1_sigma_alpha2 + 0.5 * log(2*pi) + 0.5 * log(Sigma_beta0) + 0.5 + (L/2) * log(2*pi) + (L/2) * log(det(Sigma_theta) + 10^(-10)) +
      (L/2) + (P/2) * log(2*pi) + (P/2) * log(det(Sigma_beta_confound) + 10^(-10)) +
      (P/2) + (p/2) * log(2*pi) + (p/2) * exp_ln_sigma_alpha2  + alpha_part1 +
      alpha_part2 + alpha_part3 - (n_train/2 + a_epsilon) * log(b_epsilon + 0.5 * sum((y - beta0 - X%*%beta_thresholded)^2)) +
      ((n_train/2 + a_epsilon) + 1) * exp_ln_sigma_epsilon2 +
      (b_epsilon + 0.5 * sum((y - beta0 - X%*%beta_thresholded)^2)) * exp_1_sigma_epsilon2 - (L/2 + a_beta) * log(0.5 * sum(diag(Sigma_theta) + theta^2) + b_beta) +
      (L/2 + a_beta + 1) * exp_ln_sigma_beta2 + (0.5 * sum(diag(Sigma_theta) + theta^2) + b_beta) * exp_1_sigma_beta2 -
      (p/2 + a_alpha) * log(0.5 * sum((alpha - beta)^2) + b_alpha)  + (p/2 + a_alpha + 1) * exp_ln_sigma_alpha2 +
      (0.5 * sum((alpha - beta)^2) + b_alpha) * exp_1_sigma_alpha2

    # Save ELBO and difference in ELBO values in vectors.
    ELBO_chain[iter] = ELBO
    diff = as.vector(ELBO) - as.vector(ELBO_old)
    diff_chain[iter] = diff
  
    # Save the other estimated parameters in matrices or vectors.
    beta0_chain[iter] = beta0
    beta_confound_chain[iter,] = beta_confound
    theta_chain[iter,] = theta
    beta_chain[iter,] = beta
    alpha_chain[iter,] = alpha
    beta_threshold_chain[iter,] = beta_thresholded
    alpha_threshold_chain[iter,] = alpha_thresholded
    sigma_epsilon_chain[iter] = sigma_epsilon
    sigma_alpha_chain[iter] = sigma_alpha
    sigma_beta_chain[iter] = sigma_beta

    # Save results every 10 iterations.
    if(iter %% 10 == 0){
      print(iter)
      print(diff)

      result = list()
      result$beta0 = beta0_chain
      result$theta = theta_chain
      result$beta = beta_chain
      result$beta_confound = beta_confound_chain
      result$alpha = alpha_chain
      result$sigma_epsilon = sigma_epsilon_chain
      result$sigma_beta = sigma_beta_chain
      result$sigma_alpha = sigma_alpha_chain
      result$beta_thresholded = beta_threshold_chain
      result$alpha_thresholded = alpha_threshold_chain
      result$log_likelihood = log_likelihood_chain
      result$joint = joint_chain
      result$diff_chain = diff_chain
      result$threshold = threshold_chain
      result$ELBO = ELBO_chain
      result$iter = iter
      
      save(result, file = sprintf("%sresult_rtgp_vi_left_cortex.RData", path_out))
      rm(result)

  }
  
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  time.end = Sys.time()
  time = time.end - time.start
}

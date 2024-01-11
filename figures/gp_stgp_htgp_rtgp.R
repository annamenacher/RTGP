# Script: Create figure that showcases differences between GP draw and thresholded GPs.

#################
### Libraries ###
#################

library(MASS)
library(ggplot2)
library(grid)

################
### Analysis ###
################

# Generate covariance matrix for points in `x` using given kernel function
cov_matrix <- function(x, kernel_fn, ...) {
  outer(x, x, function(a, b) kernel_fn(a, b, ...))
}

# Given x coordinates, take N draws from kernel function at those points
draw_samples <- function(x, N, seed = 1, kernel_fn, ...) {
  Y <- matrix(NA, nrow = length(x), ncol = N)
  set.seed(seed)
  for (n in 1:N) {
    K <- cov_matrix(x, kernel_fn, ...)
    Y[, n] <- mvrnorm(1, mu = rep(0, times = length(x)), Sigma = K)
  }
  Y
}

# x-coordinates
x <- seq(0, 2, length.out = 201)  
# no. of draws
N <- 1 
# for line colors
col_list <- c("black")  

# f(x)
se_kernel <- function(x, y, sigma = 1, length = 1) {
  sigma^2 * exp(- (x - y)^2 / (2 * length^2))
}
Y <- draw_samples(x, N, kernel_fn = se_kernel, length = 0.2)
df = data.frame(x = x, y = Y)
ggplot(data=df, aes(x=x, y=y)) +
  geom_line() +
  xlim(0,2) +
  ylim(-1.3, 1.3) +
  geom_hline(yintercept = 0.5, color='red') + 
  geom_hline(yintercept = -0.5, color = 'red') + 
  ggtitle("f(x)") + 
  theme(text = element_text(size = 20))       

# STGP
stgp_function = function(f, w){
  out = ifelse(abs(f) > w, sign(f)*(abs(f) - w), 0)
  return(out)
}

Y_STGP = stgp_function(Y, w = 0.5)

df = data.frame(x = x, y = Y_STGP)
ggplot(data=df, aes(x=x, y=y)) +
  geom_line() +
  xlim(0,2) +
  ylim(-1.3, 1.3) +
  geom_hline(yintercept = 0.5, color='red') + 
  geom_hline(yintercept = -0.5, color = 'red') + 
  ggtitle("STGP(f(x), 0.5)") + 
  theme(text = element_text(size = 20))       

# HTGP
htgp_function = function(f, w){
  out = ifelse(abs(f) > w, f, 0)
  return(out)
}

Y_HTGP = htgp_function(Y, w = 0.5)

df = data.frame(x = x, y1 = c(Y_HTGP[1:27,1], rep(NA, 201 - length(Y_HTGP[1:27,1]))) , y2 = c(rep(NA, 27), Y_HTGP[28:46,1], rep(NA,201-46)), y3 = c(rep(NA, 46), Y_HTGP[47:70,1], rep(NA,201-70)),
                y4 = c(rep(NA,70), Y_HTGP[71:83,1], rep(NA,201-83)), y5 = c(rep(NA,83), Y_HTGP[84:157,1], rep(NA,201-157)), y6 = c(rep(NA,157), Y_HTGP[158:169,1], rep(NA,201-169)), 
                y7 = c(rep(NA,169), Y_HTGP[170:201,1]))
ggplot(data=df, aes(x=x)) +
  geom_line(aes(y = y1)) +
  geom_line(aes(y = y2)) +
  geom_line(aes(y = y3)) +
  geom_line(aes(y = y4)) +
  geom_line(aes(y = y5)) +
  geom_line(aes(y = y6)) +
  geom_line(aes(y = y7)) +
  xlim(0,2) +
  ylim(-1.3, 1.3) +
  ylab('y') +
  geom_hline(yintercept = 0.5, color='red') + 
  geom_hline(yintercept = -0.5, color = 'red') + 
  ggtitle("HTGP(f(x), 0.5)") + 
  theme(text = element_text(size = 20))   

# RTGP
rtgp_function = function(f, w, sigma_alpha2){
  out = matrix(NA, 100,length(x))
  for(i in 1:100){
    f_tilde = rnorm(n = length(x), mean = f, sd = sqrt(sigma_alpha2))
    out[i,] = ifelse(abs(f_tilde)>0.5,1,0)*f
  }
  return(out)
}

Y_RTGP = rtgp_function(Y, w = 0.5, sigma_alpha2 = 0.01)
Y_RTGP_mean = matrix(apply(Y_RTGP, 2, mean), length(x), 1)

df = data.frame(x = x, y = Y_RTGP_mean)
ggplot(data=df, aes(x=x, y=y)) +
  geom_line() +
  xlim(0,2) +
  ylim(-1.3, 1.3) +
  geom_hline(yintercept = 0.5, color='red') + 
  geom_hline(yintercept = -0.5, color = 'red') + 
  ggtitle("RTGP(f(x), 0.5, 0.01)") + 
  theme(text = element_text(size = 20))       


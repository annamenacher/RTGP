# Script: Analyse the total variation captured with the eigendecomposition (left and right hemisphere)

#############
### Paths ###
#############

# Path to eigendecomposition (NOTE: full eigendecomposition of L = M).
path_eigen = ''

#################
### Libraries ###
#################

library(data.table)
library(ggplot2)

################
### Analysis ###
################

# Eigenvalues of eigendecomposition of left hemisphere
x = fread(paste0(path_eigen, "/eigendecomposition_L29716_phi0.17_nu1.38/eigenvalues_left_cortex.csv"), header = T)
x = x[,2]
x = as.numeric(unlist(x))

# Plot cumulative percentage of eigenvalues across number of basis functions of the eigendecomposition.
data <- data.frame(
  Eigenvalues = 1:length(x),
  Percentage = cumsum(x)/sum(x)
)

plot <- ggplot(data, aes(Eigenvalues, Percentage)) +
  xlab('Eigenvalues') +
  ylab('Cumulative Percentage') +
  theme_minimal() +
  geom_line()

print(plot)

# Sum of eigenvalues that cover at least 90% of the total variation.
which(cumsum(x)/sum(x) >= 0.9)[1]


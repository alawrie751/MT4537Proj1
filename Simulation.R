
# Description ------

# Script to complete the first part of the project (simulation section)

# Load Packages ------

library(spatstat)

# Simulation Function ------

# Function to run the simulation of a Thomas cluster process

# Inputs:
#   mu_p - mean of parent distribution (positive number)
#   mu_c - mean of child distribution (positive number)
#   sd_c - standard deviation of symmetric normal distribution of children
#          around parents (positive number)
#   rand_seed - number to set the random seed to to allow reproducibility of
#               simulations (number, has a default of 150)

# Outputs:
#   named list of 2 vectors containing n[1] and n[2] simulated data points from 
#   the specified normal distributions


ThomasSimul <- function(mu_p, mu_c, sd_c, rand_seed = 150) {
  
  # Set random seed to allow reproducibility of simulations 
  set.seed(rand_seed)
  
  # Randomly select number of parents from Poisson distribution of mean mu_p
  num_p <- rpois(1, mu_p)
  
  # Randomly generate coordinates for each of the num_p points, with equal 
  # probability everywhere in the region of [0, 10] in each direction
  x_p <- runif(num_p, 0, 10)
  y_p <- runif(num_p, 0, 10)
  #plot(x_p, y_p, xlim = c(0, 10), ylim = c(0, 10))
  
  # Convert to a spatstat ppp object
  loc_p <- ppp(x_p, y_p, window = owin(c(0,10), c(0, 10)))
  
  # Calculate the driving intensity using Gaussian densities placed at each
  # parent location with standard deviation given by sd_c
  drive_intense <- density(loc_p, sigma = sd_c, kernel = "gaussian")



}


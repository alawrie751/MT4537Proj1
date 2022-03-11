
# Description ------

# Script containing the functions to carry out the simulations for the first
# part of the project

# Load Packages ------

library(spatstat)

# Thinning Function ------

# Function to thin a simulation from a homogeneous Poisson process to an 
# inhomogeneous Poisson process simulation

# Inputs:
#   intense - intensity given by the density function from spatstat
#   full_sim - ppp object which is a simulation from a homogeneous Poisson 
#              with intensity equal to the maximum intensity seen across the
#              window
#   xylim - limit of the window in both the x and y directions (positive number)

# Outputs:
#   ppp object containing the thinned simulation from a Thomas process

thin_sim <- function(intense, full_sim, xylim) {
  
  # Create variables to store a vector containing the points being kept
  new_x <- NULL
  new_y <- NULL
  
  # Extract the x and y values into vectors and find the length of the vectors
  xvals <- full_sim$x
  yvals <- full_sim$y
  
  nvals <- length(xvals)
  
  # Convert the derived intensity function to the range [0, 1] and then convert 
  # it to a function type so that the value of the density at (x, y) is 
  # extracted easily
  prob_dens <- as.function(intense / max(intense))
  
  for (i in 1:nvals) {
    
    # Find probability density at the given point then sample a Bernoulli 
    # variable with probability of 1 equal to this probability
    prob <- prob_dens(xvals[i], yvals[i])
    ind <- rbinom(1, 1, prob)
    
    # If the indicator is 1 then keep the given point and add it to the vector
    if (ind == 1) {
      
      new_x <- c(new_x, xvals[i])
      new_y <- c(new_y, yvals[i])
      
    }
    
  }
  
  # Convert to a spatstat ppp object and return
  thinned <- ppp(new_x, new_y, window = owin(c(0, xylim), c(0, xylim)))
  
  return(thinned)
  
}


# Thomas Simulation Function ------

# Function to run the simulation of a Thomas cluster process using the 
# thinning function above

# Inputs:
#   mu_p - mean of parent distribution (positive number)
#   sd_c - standard deviation of symmetric normal distribution of children
#          around parents (positive number)
#   xylim - limit of the window in both the x and y directions (positive number,
#           has a default value of 1)
#   rand_seed - number to set the random seed to to allow reproducibility of
#               simulations (number, has an arbitrary default of 150)

# Outputs:
#   ppp object containing the simulation from a Thomas process

ThomasSimul <- function(mu_p, sd_c, xylim = 1, rand_seed = 150) {
  
  # Error traps to ensure that user inputs are valid
  if (!is.numeric(mu_p)) stop("invalid arguments. mu_p must be a numeric")
  if (!is.numeric(sd_c)) stop("invalid arguments. sd_c must be a numeric")
  if (!is.numeric(xylim)) stop("invalid arguments. xylim must be a numeric")
  if (!is.numeric(rand_seed)) stop("invalid arguments. rand_seed must be a numeric")
  
  if (mu_p <= 0) stop("invalid arguments. mu_p must be > 0")
  if (sd_c <= 0) stop("invalid arguments. sd_c must be > 0")
  if (xylim <= 0) stop("invalid arguments. xylim must be > 0")
  
  # Set random seed to allow reproducibility of simulations 
  set.seed(rand_seed)
  
  # Randomly select number of parents from Poisson distribution of mean mu_p
  num_p <- rpois(1, mu_p)
  
  # Randomly generate coordinates for each of the num_p points, with equal 
  # probability everywhere in the region of [0, xylim] in each direction
  x_p <- runif(num_p, 0, xylim)
  y_p <- runif(num_p, 0, xylim)
  #plot(x_p, y_p, xlim = c(0, xylim), ylim = c(0, xylim))
  
  # Convert to a spatstat ppp object
  loc_p <- ppp(x_p, y_p, window = owin(c(0, xylim), c(0, xylim)))
  
  # Calculate the driving intensity using Gaussian densities placed at each
  # parent location with standard deviation given by sd_c
  drive_intense <- density(loc_p, sigma = sd_c, kernel = "gaussian")
  
  # Create sample of a homogeneous Poisson process with density equal to the
  # highest seen in the derived driving intensity above
  large_sim <- rpoispp(max(drive_intense), win = owin(c(0, xylim), c(0, xylim)))
  
  # Thin this using the thinning function described above
  thinned_sim <- thin_sim(drive_intense, large_sim, xylim)
  
  return(list(thin_sim = thinned_sim, intense = drive_intense))
  
}


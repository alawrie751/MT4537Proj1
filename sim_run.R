
# Description ------

# Script to run the simulations for the first part of the project, create the 
# plots and fit the models as required

# Source the Functions ------

# Note that this also loads the spatstat library
source("sim_funcs.R")

# Create 2 Simulations, Clip to Window and Fit Models ------

sim1 <- ThomasSimul(30, 0.1, 200)
sim2 <- ThomasSimul(45, 0.05, 465)


## Clip to window ##

thom_sim1 <- kppm(sim1, clusters = "Thomas", method = "mincon")
mat_sim1 <- kppm(sim1, clusters = "MatClust", method = "mincon")

thom_sim2 <- kppm(sim2, clusters = "Thomas", method = "mincon")
mat_sim2 <- kppm(sim2, clusters = "MatClust", method = "mincon")




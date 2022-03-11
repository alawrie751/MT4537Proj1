
# Description ------

# Script to run the simulations for the first part of the project, create the 
# plots and fit the models as required

# Source the Functions ------

# Note that this also loads the spatstat library
source("src/sim_funcs.R")

# Create 2 Simulations, Clip to Window and Fit Models ------

sim1 <- ThomasSimul(30, 0.1, rand_seed = 200)
sim2 <- ThomasSimul(45, 0.05, rand_seed = 465)

# Save plots of these simulations and the driving intensity
pdf(file = "results/simplots/plot_sim1.pdf")
plot(sim1$thin_sim, main = "Plot of the First Simulation from a 
     Thomas Cluster Process")
dev.off()

pdf(file = "results/simplots/intense_sim1.pdf")
plot(sim1$intense, main = "Plot of the Driving Intensity for the 
     First Simulation from a Thomas Cluster Process")
dev.off()

pdf(file = "results/simplots/plot_sim2.pdf")
plot(sim2$thin_sim, main = "Plot of the Second Simulation from a 
     Thomas Cluster Process")
dev.off()

pdf(file = "results/simplots/intense_sim2.pdf")
plot(sim2$intense, main = "Plot of the Driving Intensity for the 
     Second Simulation from a Thomas Cluster Process")
dev.off()

## Clip to window ##

thom_sim1 <- kppm(sim1, clusters = "Thomas", method = "mincon")
mat_sim1 <- kppm(sim1, clusters = "MatClust", method = "mincon")

thom_sim2 <- kppm(sim2, clusters = "Thomas", method = "mincon")
mat_sim2 <- kppm(sim2, clusters = "MatClust", method = "mincon")





# Description ------

# Script to run the simulations for the first part of the project, 
# create the plots and fit the models as required

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

# Clip to Nominal Window ------

# For simulation 1 clip to [0.2, 0.8] on each axis

# Get vector of indicies to remove
indremove <- NULL

for (i in 1:length(sim1$thin_sim$x)) {
  
  xi <- sim1$thin_sim$x[i]
  yi <- sim1$thin_sim$y[i]
  
  if (xi < 0.2 | xi > 0.8) {
    
    indremove <- c(indremove, i) 
    
  } else if (yi < 0.2 | yi > 0.8) {
    
    indremove <- c(indremove, i)
    
  }
  
}

# Remove from ppp object
sim1$thin_sim$x <- sim1$thin_sim$x[-indremove]
sim1$thin_sim$y <- sim1$thin_sim$y[-indremove]

pdf(file = "results/simplots/clip_sim1.pdf")
plot(sim1$thin_sim, main = "Plot of the First Simulation Clipped")
dev.off()

# For simulation 2 clip to [0.10, 0.90] on each axis

# Get vector of indicies to remove
indremove2 <- NULL

for (i in 1:length(sim2$thin_sim$x)) {
  
  xi <- sim2$thin_sim$x[i]
  yi <- sim2$thin_sim$y[i]
  
  if (xi < 0.1 | xi > 0.9) {
    
    indremove2 <- c(indremove2, i) 
    
  } else if (yi < 0.1 | yi > 0.9) {
    
    indremove2 <- c(indremove2, i)
    
  }
  
}

# Remove from ppp object
sim2$thin_sim$x <- sim2$thin_sim$x[-indremove2]
sim2$thin_sim$y <- sim2$thin_sim$y[-indremove2]

pdf(file = "results/simplots/clip_sim2.pdf")
plot(sim2$thin_sim, main = "Plot of the Second Simulation Clipped")
dev.off()

# Fit models ------

thom_sim1 <- kppm(sim1$thin_sim, 
                  clusters = "Thomas", 
                  method = "mincon")
mat_sim1 <- kppm(sim1$thin_sim, 
                 clusters = "MatClust", 
                 method = "mincon")

thom_sim2 <- kppm(sim2$thin_sim, 
                  clusters = "Thomas", 
                  method = "mincon")
mat_sim2 <- kppm(sim2$thin_sim, 
                 clusters = "MatClust", 
                 method = "mincon")


datPCF_tom1 <- pcfinhom(sim1$thin_sim, 
                        thom_sim1, 
                        correction = "translate")
datPCF_sim_tom1 <- pcfinhom(simulate(thom_sim1, nsim = 1)[[1]], 
                            thom_sim1, 
                            correction = "translate")

datPCF_mat1 <- pcfinhom(sim1$thin_sim, 
                        mat_sim1, 
                        correction = "translate")
datPCF_sim_mat1 <- pcfinhom(simulate(mat_sim1, nsim = 1)[[1]], 
                            mat_sim1, 
                            correction = "translate")

datPCF_tom2 <- pcfinhom(sim2$thin_sim, 
                        thom_sim2, 
                        correction = "translate")
datPCF_sim_tom2 <- pcfinhom(simulate(thom_sim2, nsim = 1)[[1]], 
                            thom_sim2, 
                            correction = "translate")

datPCF_mat2 <- pcfinhom(sim1$thin_sim, 
                        mat_sim2, 
                        correction = "translate")
datPCF_sim_mat2 <- pcfinhom(simulate(mat_sim2, nsim = 1)[[1]], 
                            mat_sim2, 
                            correction = "translate")

pdf(file = "results/simplots/all_pcfs.pdf")
par(mfrow = c(2, 3))
plot(datPCF_tom1, main = "Real Data (1)")
plot(datPCF_sim_tom1, main = "Thomas (1)")
plot(datPCF_sim_mat1, main = "Matern (1)")
plot(datPCF_tom2, main = "Real Data (2)")
plot(datPCF_sim_tom2, main = "Thomas (2)")
plot(datPCF_sim_mat2, main = "Matern (2)")
dev.off()
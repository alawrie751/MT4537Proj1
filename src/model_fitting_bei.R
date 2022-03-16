
# Description ------

# Script to run the modelling for the bei data in the 
# second part of the project

# Load packages and bei Data------

library(spatstat)

dat <- bei
dat_extra <- bei.extra

pdf(file = "results/initialintenseplots/beidat.pdf")
plot(bei, main = "bei data")
dev.off()

# set the seed so simulations are reproducible
set.seed(150)

# Fit Kernel Smooth ------

# Fit kernel smooth as initial estimate of intensity
initial_est <- density(dat)

# Save plot of this intensity
pdf(file = "results/initialintenseplots/initial_density.pdf")
plot(initial_est, 
     main = "Initial intensity Estimate Using Function Defaults")
dev.off()

# Default kernel is Gaussian so this is use since it 
# hasn't been overridden

# the bandwidth sigma is found by a "simple rule of thumb" 
# relating only to the size of the window (assuming some 
# 0.1 * shortest side length)

# Fit densities with a large and small bandwidth to give 
# markedly different estimates of the intensity
large_band <- density(dat, sigma = 200)
small_band <- density(dat, sigma = 20)

pdf(file = "results/initialintenseplots/initial_large.pdf")
plot(large_band)
dev.off()

pdf(file = "results/initialintenseplots/initial_small.pdf")
plot(small_band)
dev.off()

# Fit Models Using ppm for x and y ------

mod_stat <- ppm(dat ~ 1)

mod_x_y <- ppm(dat ~ x + y)

mod_xy <- ppm(dat ~ x * y)

mod_poly2_xy <- ppm(dat ~ poly(x, 2) + poly(y, 2))

mod_poly3_xy <- ppm(dat ~ poly(x, 3) + poly(y, 3))

mod_poly6_xy <- ppm(dat ~ poly(x, 6) + poly(y, 6))

par(mfrow = c(2, 2), mar = c(0.5, 0.5, 0.5, 0.5))
plot(predict(mod_poly2_xy))
plot(predict(mod_poly3_xy))
plot(predict(mod_poly6_xy))
plot(density(dat, sigma = 50))

pdf(file = "results/ppmmodels/best_xy.pdf")
plot(predict(mod_poly6_xy), 
     main = "Model Fitted with Polynomial of Degree 6")
dev.off()

datPCF_x <- pcfinhom(dat, 
                     mod_poly6_xy, 
                     correction = "translate")
datPCF_sim_x <- pcfinhom(simulate(mod_poly6_xy, nsim = 1)[[1]],
                         mod_poly6_xy, 
                         correction = "translate")

pdf(file = "results/ppmmodels/pcf_realxy.pdf")
plot(datPCF_x, main = "Real Data")
dev.off()

pdf(file = "results/ppmmodels/pcf_simxy.pdf")
plot(datPCF_sim_x, main = "Simulated Data")
dev.off()

# Fit Models Using ppm for Elevation and Gradient ------

mod_elev <- ppm(dat ~ elev, data = dat_extra)

mod_grad <- ppm(dat ~ grad, data = dat_extra)

mod_g_e <- ppm(dat ~ elev + grad, data = dat_extra)

mod_ge <- ppm(dat ~ elev * grad, data = dat_extra)

pdf(file = "results/ppmmodels/best_ge.pdf")
plot(predict(mod_ge), 
     main = "Model Fitted with the Interaction Term")
dev.off()

datPCF <- pcfinhom(dat, 
                   mod_ge, 
                   correction = "translate")
datPCF_sim <- pcfinhom(simulate(mod_ge, nsim = 1)[[1]], 
                       mod_ge, 
                       correction = "translate")

pdf(file = "results/ppmmodels/pcf_real.pdf")
plot(datPCF, main = "Real Data")
dev.off()

pdf(file = "results/ppmmodels/pcf_sim.pdf")
plot(datPCF_sim, main = "Simulated Data")
dev.off()

# Fit Thomas Models ------

mod_tom <- kppm(dat ~ poly(x, 6) + poly(y, 6), 
                clusters = "Thomas", 
                method = "palm")

pdf(file = "results/ppmmodels/best_tom.pdf")
plot(predict(mod_tom), 
     main = "Thomas Cluster Model Fitted with Polynomial of Degree 6")
dev.off()

pdf(file = "results/ppmmodels/best_tom_dat.pdf")
plot(mod_tom, 
     main = "Thomas Cluster Model Fitted with Polynomial of Degree 6 with Data")
dev.off()

datPCF_tom <- pcfinhom(dat, 
                       mod_tom, 
                       correction = "translate")
datPCF_sim_tom <- pcfinhom(simulate(mod_tom, nsim = 1)[[1]], 
                           mod_tom, 
                           correction = "translate")

pdf(file = "results/ppmmodels/pcf_realtom.pdf")
plot(datPCF_tom, main = "Real Data")
dev.off()

pdf(file = "results/ppmmodels/pcf_simtom.pdf")
plot(datPCF_sim_tom, main = "Simulated Data")
dev.off()

# Simulations from the Cox Model ------

# Set a new seed for these simulations
set.seed(500)

# Simulate from the Thomas model fitted above
tomsim <- simulate(mod_tom, nsim = 2)

pdf(file = "results/ppmmodels/tomsim.pdf")
plot(tomsim, 
     main = "Data Simulated from the Thomas Model")
dev.off()

# Simulations from the Poisson Model ------

set.seed(600)

poly6_sim <- simulate(mod_poly6_xy, nsim = 2)

pdf(file = "results/ppmmodels/poly6sim.pdf")
plot(poly6_sim, 
     main = "Data Simulated from the Poisson Process Model")
dev.off()

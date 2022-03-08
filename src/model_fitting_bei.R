
# Description ------

# Script to run the modelling for the bei data in the second part of the project

# Load packages and bei Data------

library(spatstat)

dat <- bei
dat_extra <- bei.extra

# Fit Kernel Smooth ------

# Fit kernel smooth as initial estimate of intensity
initial_est <- density(dat)

# Save plot of this intensity
pdf(file = "results/initialintenseplots/initial_density.pdf")
plot(initial_est, main = "Initial intensity Estimate Using Function Defaults")
dev.off()

# Default kernel is gaussian so this is use since it hasn't been overridden

# the bandwidth sigma is found by a "simple rule of thumb" relating only to
# the size of the window (assuming some 0.1 * shortest side length)

# Fit densities with a large and small bandwidth to give markedly different
# estimates of the intensity
large_band <- density(dat, sigma = 200)
small_band <- density(dat, sigma = 20)

pdf(file = "results/initialintenseplots/initial_large.pdf")
plot(large_band)
dev.off()

pdf(file = "results/initialintenseplots/initial_small.pdf")
plot(small_band)
dev.off()








# Code used for the example in Section 1.1.2 and 1.1.3

library(spatstat)

set.seed(123)

num_p <- rpois(1, 11)

x_p <- runif(num_p, 0, 1)
y_p <- runif(num_p, 0, 1)

pdf(file = "results/simplots/par_eg1.pdf")
plot(x_p, y_p, xlim = c(0, 1), ylim = c(0, 1), xlab = "x", ylab = "y", 
     main = "Plot of Parent Points for the Example")
dev.off()

loc_par <- ppp(x_p, y_p, window = owin(c(0, 1), c(0, 1)))

drive_intense <- density(loc_par, sigma = 0.1, kernel = "gaussian")

pdf(file = "results/simplots/intense_eg1.pdf")
plot(drive_intense, main = "Plot of the Driving Intensity for the 
     Example")
dev.off()

full_sim <- rpoispp(max(drive_intense), win = owin(c(0, 1), c(0, 1)))

pdf(file = "results/simplots/fullsim_eg2.pdf")
plot(full_sim, main = "Plot of the Full Homogeneous Poisson Simulation")
dev.off()

new_x <- NULL
new_y <- NULL

xvals <- full_sim$x
yvals <- full_sim$y

nvals <- length(xvals)

prob_dens <- as.function(drive_intense / max(drive_intense))

for (i in 1:nvals) {
  
  prob <- prob_dens(xvals[i], yvals[i])
  ind <- rbinom(1, 1, prob)
  
  if (ind == 1) {
    
    new_x <- c(new_x, xvals[i])
    new_y <- c(new_y, yvals[i])
    
  }
  
}

thinned <- ppp(new_x, new_y, window = owin(c(0, 1), c(0, 1)))

pdf(file = "results/simplots/thinsim_eg2.pdf")
plot(full_sim, main = "Plot of the Full Homogeneous Poisson Simulation")
points(thinned, pch = 16, col = "red")
dev.off()

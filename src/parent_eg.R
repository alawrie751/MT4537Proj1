
# Code used for the example in Section 1.1.2

library(statspat)

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

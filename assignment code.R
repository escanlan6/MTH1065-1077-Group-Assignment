# M = number of paths
# N =  T /h (number of steps, where h is the step size)

# Write a function SVMil which takes as input M (as above) and N (as above), and outputs
# an (N + 1) × M matrix in which each column is a simulated path of the stock price process
# {St : t ≥ 0} in the market described above using the Milstein scheme. Explain and justify
# any corrections or adjustments you make to the standard Milstein scheme.
N <- 10
M <- 5

S_0 <- 110
V_0 <- 0.3
r <- 0.04
xi <- 0.5
theta <- 0.06
kappa <- 1.8
rho <- -0.2
maturity <- 2
h <- maturity / N

# Brownian motion function
BM <- function(rho, h, M = 1) {
  # Generate 2 vectors of standard normal random variables
  Z1 <- rnorm(M)
  Z2 <- rnorm(M)
  # Correlated Brownian increments
  W_v <- sqrt(h) * Z1
  W_s <- sqrt(h) * (rho * Z1 + sqrt(1 - rho^2) * Z2)
  # Return brownian motions
  list(W_v = W_v, W_s = W_s)
}

# Milstein function
SVMil <- function(M, N) {
  #stuff
}

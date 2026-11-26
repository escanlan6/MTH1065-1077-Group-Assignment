# Defining all variables
S0 <- 110
V0 <- 0.3
r <- 0.04
xi <- 0.5
theta <- 0.06
kappa <- 0.8
ro <- -0.2

BU <- 130

maturity <- 2

# Defining SVEM function to calculate V_t and then S_t within itself.

SVEM <- function(S0, V0, r, xi, theta, kappa, ro, M, N) {
  h <- maturity/N
  t <- seq(0, maturity, by = h)

  #Create S and V to calculate concurrently, since the V_t value is needed to
  #calculate the S_t value
  S_mat <- matrix(NA, nrow = N+1, ncol = M)
  V_mat <- matrix(NA, nrow = N+1, ncol = M)

  S_mat[1, ] <- S0
  V_mat[1, ] <- V0

  #Correlated Brownian motions
  dWv <- matrix(rnorm(N * M, mean = 0, sd = sqrt(h)), nrow = N, ncol = M)
  dZ <- matrix(rnorm(N * M, mean = 0, sd = sqrt(h)), nrow = N, ncol = M)
  dWs <- ro * dWv + sqrt(1 - ro^2)*dZ

  for (i in 1:N) {
    
    V_i <- V_mat[i, ]
    S_i <- S_mat[i, ]

    V_mat[i+1, ] <- V_i + kappa * (theta - V_i) * h + xi * (V_i^2) * dWv[i, ]

    S_mat[i+1, ] <- S_i + r * S_i * h + sqrt(V_i) * S_i * dWs[i, ]

  }
  return(S_mat)

}

test_N_a_i <- 20
test_M_a_i <- 20
test_sample_a_i <- SVEM(S0, V0, r, xi, theta, kappa, ro, M, N)

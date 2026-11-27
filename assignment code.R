# Simulation for finance assignment


# Defining all variables
S0 <- 110
V0 <- 0.3
r <- 0.04
xi <- 0.5
theta <- 0.06
kappa <- 0.8
rho <- -0.2

K <- 110

BU <- 130

maturity <- 2

# Brownian motion function
BM <- function(rho, maturity, N, M) {
  # Define step size
  h <- maturity / N
  # Generate 2 vectors of standard normal random variables
  Z1 <- rnorm(N*M)
  Z2 <- rnorm(N*M)
  # Correlated Brownian motions
  W_v <- sqrt(h) * Z1
  W_s <- sqrt(h) * (rho * Z1 + sqrt(1 - rho^2) * Z2)
  # Return brownian motions
  list(W_v = W_v, W_s = W_s)

}

# Euler-Maruyama scheme
# Defining SVEM function to calculate V_t and then S_t within itself.
SVEM <- function(S0, V0, r, xi, theta, kappa, rho, M, N) {
  h <- maturity / N
  t <- seq(0, maturity, by = h)

  # Create S and V to calculate concurrently, since the V_t value is needed to
  # Calculate the S_t value
  S_mat <- matrix(NA, nrow = N+1, ncol = M)
  V_mat <- matrix(NA, nrow = N+1, ncol = M)

  S_mat[1, ] <- S0
  V_mat[1, ] <- V0

  # Correlated Brownian motions
  dW <- BM(rho, maturity, N, M)
  dWv <- matrix(dW$W_v, nrow = N, ncol = M)
  dWs <- matrix(dW$W_s, nrow = N, ncol = M)

  for (i in 1:N) {
    
    V_i <- V_mat[i, ]
    S_i <- S_mat[i, ]

    # V_pos is V_i if V_i is non-negative, and is 0 if negative
    V_pos <- max(V_i, 0)

    V_mat[i+1, ] <- V_i + kappa * (theta - V_i) * h + xi * (V_i^2) * dWv[i, ]

    S_mat[i+1, ] <- S_i + r * S_i * h + sqrt(V_pos) * S_i * dWs[i, ]

  }
  return(S_mat)

}

test_N_a_i <- 20
test_M_a_i <- 20
test_sample_a_i <- SVEM(S0, V0, r, xi, theta, kappa, rho, test_M_a_i, test_N_a_i)

# Adjustments to Euler-Maruyama function: since sqrt(V_t) is not defined for V_t<0
# For the calculation of S_t, sqrt(max(V_t, 0)) is used instead.

# Milstein scheme function
SVMil <- function(S0, V0, r, xi, theta, kappa, rho, M, N) {
  h <- maturity / N
  t <- seq(0, maturity, by = h)

  #Create S and V to calculate concurrently, since the V_t value is needed to
  #calculate the S_t value
  S_mat <- matrix(NA, nrow = N+1, ncol = M)
  V_mat <- matrix(NA, nrow = N+1, ncol = M)

  S_mat[1, ] <- S0
  V_mat[1, ] <- V0

  # Correlated Brownian motions
  dW <- BM(rho, maturity, N, M)
  dWv <- matrix(dW$W_v, nrow = N, ncol = M)
  dWs <- matrix(dW$W_s, nrow = N, ncol = M)

  for (i in 1:N) {
    
    V_i <- V_mat[i, ]
    S_i <- S_mat[i, ]

    V_pos <- max(V_i, 0)

    # Milstein addition to V
    V_mat[i+1, ] <- V_i + kappa * (theta - V_i) * h + xi * (V_i^2) * dWv[i, ] + xi^2 * (V_i^3) * (dWv[i, ]^2 - h)

    # Milstein addition to S
    S_mat[i+1, ] <- S_i + r * S_i * h + sqrt(V_pos) * S_i * dWs[i, ] + 0.5 * V_i * S_i * (dWs[i, ]^2 - h)

  }
  return(S_mat)

}


# (b) (i)

# Sample size sets for N and M
N_list <- c(5, 10, 100, 1000)
len_n <- length(N_list)
M_list <- c(10, 100, 1000, 5000, 10000)
len_m <- length(M_list)

mc_mat_EM <- matrix(data=NA, nrow = len_n, ncol = len_m)

mc_mat_Mil <- matrix(data=NA, nrow = len_n, ncol = len_m)


for (i in 1:len_n){
  for (j in 1:len_m) {
    EM <- SVEM(S0, V0, r, xi, theta, kappa, rho, M_list[j], N_list[i])
    Mil <- SVMil(S0, V0, r, xi, theta, kappa, rho, M_list[j], N_list[i])

    S_T_EM <- EM[nrow(EM), ]
    S_T_Mil <- Mil[nrow(Mil), ]

    p_EM <- pmax(K - S_T_EM, 0)
    p_Mil <- pmax(K - S_T_Mil, 0)
    
    # Discounted payoff:
    mc_mat_EM[i, j] <- exp(-r * maturity) * mean(p_EM)
    mc_mat_Mil[i, j] <- exp(-r * maturity) * mean(p_Mil)

  }
}

# Making data frame of results for cleaner tabulation
mc_results_EM <- as.data.frame(mc_mat_EM)
rownames(mc_results_EM) <- c("5 steps", "10 steps", "100 steps", "1000 steps")
colnames(mc_results_EM) <- c("10 paths", "100 paths", "1000 paths", "5000 paths", "10000 paths")

mc_results_Mil <- as.data.frame(mc_mat_Mil)
rownames(mc_results_Mil) <- c("5 steps", "10 steps", "100 steps", "1000 steps")
colnames(mc_results_Mil) <- c("10 paths", "100 paths", "1000 paths", "5000 paths", "10000 paths")

print(mc_results_EM)
print(mc_results_Mil)

# There is a higher prevalence of results = 0 for the 10 paths, particularly in the
# Euler-Maruyama scheme

# There is overall agreement of results from 21-27 in both methods at the larger sample sizes.


# (b) (ii)

# Black-Scholes version calculation (sigma = sqrt(V_i) = sqrt(theta))
d1 = (log(S0/K) + maturity * (r + theta / 2)) / (sqrt(theta) * sqrt(maturity))
d2 = d1 - sqrt(theta) * sqrt(maturity)
# European Put option in the Black-Scholes market
put = pnorm(-d2) * K *exp(-r * maturity) - pnorm(-d1) * S0
put


#(c) 

M_c <- 10000
N_c <- 1000

mat_c_EM <- SVEM(S0, V0, r, xi, theta, kappa, rho, M_c, N_c)
mat_c_Mil <- SVMil(S0, V0, r, xi, theta, kappa, rho, M_c, N_c)

payoff_mat <- function(S_t, K, BU, r, maturity) {
  M_T <- apply(S_t, 2, max)
  S_T <- S_t[nrow(S_t), ]
  if (M_T < BU) {
    pay <- exp(-r * maturity) * pmax(K - S_t, 0)
  }
  return(pay)
}

payoff_c_EM <- payoff_mat(mat_c_EM)
payoff_c_Mil <- payoff_mat(mat_c_Mil)

avg_price_c_EM <- mean(payoff_c_EM)
avg_price_c_Mil <- mean(payoff_c_Mil)

# # (d)

# Sample size sets for N and M
N <- 1000
M0 <- 5000

EM <- SVEM(S0, V0, r, xi, theta, kappa, rho, M0, N)
Mil <- SVMil(S0, V0, r, xi, theta, kappa, rho, M0, N)

S_T_EM <- EM[nrow(EM), ]
S_T_Mil <- Mil[nrow(Mil), ]

# Discounted payoffs
Y_EM <- payoff_mat(EM)
Y_Mil <- payoff_mat(Mil)

# Estimated optimal control coefficient: b* = Cov(Y, S_T) / Var(S_T)
b_star_EM <- cov(Y_EM, S_T_EM) / var(S_T_EM)
b_star_Mil <- cov(Y_Mil, S_T_Mil) / var(S_T_Mil)

b_star_EM
b_star_Mil


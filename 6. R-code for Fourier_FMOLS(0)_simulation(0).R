#++++++++++++++++++++++++++++++++++++++++++++++++
# MASTER SCRIPT: Fourier-FMOLS Monte Carlo (1000 Reps)
#++++++++++++++++++++++++++++++++++++++++++++++++
setwd("C:/R/fmols")
library(cointReg)

# --- Simulation Parameters ---
T_val <- 200     # Sample size
reps <- 1000     # Number of repetitions
true_beta <- 2.5 
bandwidth <- 8   # Fixed bandwidth for comparison

# Matrices to store t-statistics and slopes for 4 scenarios
sim_slopes <- matrix(NA, nrow = reps, ncol = 4)
sim_tstats <- matrix(NA, nrow = reps, ncol = 4)
colnames(sim_slopes) <- colnames(sim_tstats) <- c("Strong", "Weak", "Very_Weak", "None")

#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 1: Refined Fourier FMLOS Function
#++++++++++++++++++++++++++++++++++++++++++++++++
fourier_fmols_manual <- function(y, x, bandwidth, k_candidates = 1:5) {
  y <- as.matrix(y)
  x0 <- as.matrix(x)
  T_orig <- nrow(y)
  t_index <- 1:T_orig

  best_aic <- Inf
  best_k <- NA
  best_ols_u <- NULL
  best_X_full <- NULL

  # --- AIC Selection for k ---
  for (k in k_candidates) {
    sin_term <- sin(2 * pi * k * t_index / T_orig)
    cos_term <- cos(2 * pi * k * t_index / T_orig)
    X_full_k <- cbind(1, sin_term, cos_term, x0)

    beta_ols <- qr.solve(X_full_k, y)
    u_hat <- y - X_full_k %*% beta_ols
    ssr <- sum(u_hat^2)
    current_aic <- T_orig * log(ssr / T_orig) + 2 * ncol(X_full_k)

    if (current_aic < best_aic) {
      best_aic <- current_aic
      best_k <- k
      best_ols_u <- u_hat
      best_X_full <- X_full_k
    }
  }

  # --- FMOLS Math ---
  dx <- diff(x0)
  u_hat_trim <- best_ols_u[-1, , drop = FALSE]
  x_trim <- best_X_full[-1, , drop = FALSE]
  y_trim <- y[-1, , drop = FALSE]
  T_trim <- nrow(x_trim)

  w <- cbind(u_hat_trim, dx)
  Sigma <- (t(w) %*% w) / T_trim
  Lambda <- matrix(0, 2, 2)

  for (j in 1:bandwidth) {
    weight <- 1 - (j / bandwidth)
    Gamma_j <- (t(w[(j+1):T_trim, ]) %*% w[1:(T_trim-j), ]) / T_trim
    Lambda <- Lambda + weight * Gamma_j
  }

  Omega <- Sigma + Lambda + t(Lambda)
  Delta <- Sigma + Lambda

  # Corrections
  y_plus <- y_trim - dx %*% solve(Omega[2,2]) %*% Omega[1,2]
  Delta_plus <- Delta[1,2] - Omega[1,2] %*% solve(Omega[2,2]) %*% Delta[2,2]
  correction <- matrix(c(0, 0, 0, T_trim * as.numeric(Delta_plus)), 4, 1)

  # Estimate Theta
  cross_X <- t(x_trim) %*% x_trim
  theta <- qr.solve(cross_X, t(x_trim) %*% y_plus - correction)
  
  # Standard Error Calculation
  omega_0_x <- Omega[1,1] - Omega[1,2] %*% solve(Omega[2,2]) %*% Omega[2,1]
  var_covar <- as.numeric(omega_0_x) * solve(cross_X)
  se_slope <- sqrt(var_covar[4, 4])
  t_stat <- (theta[4, 1] - true_beta) / se_slope # Testing if Beta = 2.5

  return(list(slope = theta[4, 1], tstat = t_stat))
}

#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 2: The Simulation Loop
#++++++++++++++++++++++++++++++++++++++++++++++++
set.seed(42)
cat("Starting 1000 Repetitions...\n")

for (i in 1:reps) {
  # Generate X (Random Walk)
  x <- cumsum(rnorm(T_val))
  
  # Generate 4 Scenarios
  y1 <- 10 + true_beta * x + rnorm(T_val)                            # Strong
  y2 <- 10 + true_beta * x + arima.sim(list(ar=0.8), T_val)          # Weak
  y3 <- 10 + true_beta * x + arima.sim(list(ar=0.95), T_val)         # Very Weak
  y4 <- 10 + true_beta * x + cumsum(rnorm(T_val))                    # None
  
  for (j in 1:4) {
    y_current <- get(paste0("y", j))
    res <- fourier_fmols_manual(y_current, x, bandwidth)
    sim_slopes[i, j] <- res$slope
    sim_tstats[i, j] <- res$tstat
  }
  if (i %% 100 == 0) cat("Repetition:", i, "\n")
}

#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 3: Analysis of Results
#++++++++++++++++++++++++++++++++++++++++++++++++
# Rejection Rate at 5% (Testing H0: Beta = 2.5)
# If H0 is TRUE, this should be 0.05.
power_size <- apply(sim_tstats, 2, function(x) mean(abs(x) > 1.96))

cat("\nEmpirical Rejection Rates (H0: Beta = 2.5):\n")
print(power_size)

# Summary of Slopes
cat("\nSummary of Slope Estimates:\n")
print(round(apply(sim_slopes, 2, summary), 4))
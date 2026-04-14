#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 1: set a working directory and library
#++++++++++++++++++++++++++++++++++++++++++++++++
setwd("C:/R/fmols")
library(cointReg)

T_val <- 200
mout <- matrix(data=NA, nrow=4, ncol=1)

#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 2: load data
#++++++++++++++++++++++++++++++++++++++++++++++++
my_data <- read.csv("data(4).csv")

for (j in 1:4) {
  # Standardize indexing to 2:T (199 observations)
  if (j == 1) {
    y <- as.matrix(my_data[2:T_val, 3])
    x <- as.matrix(my_data[2:T_val, 2])
    k_val <- 4
  } else if (j == 2) {
    y <- as.matrix(my_data[2:T_val, 4])
    x <- as.matrix(my_data[2:T_val, 2])
    k_val <- 2
  } else if (j == 3) {
    y <- as.matrix(my_data[2:T_val, 5])
    x <- as.matrix(my_data[2:T_val, 2])
    k_val <- 1
  } else {
    y <- as.matrix(my_data[2:T_val, 6])
    x <- as.matrix(my_data[2:T_val, 2])
    k_val <- 1
  }

  #++++++++++++++++++++++++++++++++++++++++++++++++
  # Step 3: FMOLS analysis
  #++++++++++++++++++++++++++++++++++++++++++++++++
  T_total <- nrow(y)
  t_index <- 1:T_total

  # Fourier Terms
  sin_term <- sin(2 * pi * k_val * t_index / T_total)
  cos_term <- cos(2 * pi * k_val * t_index / T_total)

  # FIX 1: Changed 'cbine' to 'cbind'
  # FIX 2: Added 'intercept' name to the column of 1s
  deter_mat <- cbind(intercept = 1, sin_term, cos_term)

  # FMOLS Estimation
  fmols_model <- cointRegFM(y = y, x = x, deter = deter_mat, kernel = "ba", bandwidth = 5)

  # FIX 3: theta[4] is the X coefficient
  # (theta[1]=Int, [2]=Sin, [3]=Cos, [4]=X)
  fmols_coeff <- fmols_model$theta[4]

  mout[j,1] <- fmols_coeff

  cat("\n--- Scenario", j, " (k =", k_val, ") ---\n")
  print(fmols_model)
}

cat("\nFinal Results (X-Slopes):\n")
print(mout)

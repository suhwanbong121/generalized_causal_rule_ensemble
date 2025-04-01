library(SuperLearner)
library(splines)
library(grpreg)
library(tsbart)
library(bart)
library(dbarts)

# Logistic function




n <- 500  # Number of individuals for simulation


# Data generating function, conditional on X




apply(fin, 2, mean)


true_estimator <- function(z_grid, X, Z, Y) {
  n <- nrow(X)
  L <- length(z_grid)
  
  f_x <- X[, 1] + X[, 3] + X[, 4]
  r1_x <- as.numeric(X[, 1] == 1 & X[, 2] == 0)
  r2_x <- as.numeric(X[, 5] == 1 & X[, 6] == 0)
  
  # Repeat f_x, r1_x, r2_x across columns (one per z)
  f_mat <- matrix(f_x, nrow = n, ncol = L)
  r1_mat <- matrix(r1_x, nrow = n, ncol = L)
  r2_mat <- matrix(r2_x, nrow = n, ncol = L)
  
  # Matrix of z values: each column is one z
  z_mat <- matrix(rep(z_grid, each = n), nrow = n, ncol = L)
  
  # Compute true mu_i(z) for each i and z
  mu_true_matrix <- f_mat + z_mat * r1_mat - z_mat * r2_mat
  
  return(mu_true_matrix)
}

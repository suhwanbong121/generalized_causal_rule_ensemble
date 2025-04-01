#### Data generation
expit <- function(x) {
  1 / (1 + exp(-x))
}

# Propensity function
pi_x <- function(x) {
  expit(1 - x[, 1] + x[, 4] - x[, 5])
}

# Simulate covariates
simulate_X <- function(n = 1000) {
  mat <- matrix(rbinom(n * 6, size = 1, prob = 0.5), nrow = n, ncol = 6)
  colnames(mat) <- paste0("X", 1:6)
  mat
}

generate_sim_data <- function(X, phi = 1, sigma = 0.5, k = 3) {
  n <- nrow(X)
  pi_vals <- pi_x(X)
  Z <- rbeta(n, shape1 = phi * pi_vals, shape2 = phi * (1 - pi_vals))
  
  f_x <- X[, 2] + X[, 3] + X[, 6]
  f_z <- 1 + Z + Z^2 - Z^3
  r1_x <- as.numeric(X[, 1] == 1 & X[, 2] == 0)
  r2_x <- as.numeric(X[, 5] == 1 & X[, 6] == 0)
  
  mu <- f_z + 2 * k * (Z * r1_x - Z * r2_x)
  Y <- rnorm(n, mean = mu, sd = sigma)
  
  data.frame(X, Z = Z, Y = Y)
}

# ---- Discovery Step: tsBART-based Rule Selection ---- #
rule_generation <- function(data, z_grid = NULL, J = 10, lambda = 0.1, h = 0.05, nburn = 50, nsim = 10, ntree = 40, method) {
  library(SuperLearner)
  library(splines)
  library(grpreg)
  library(tsbart)
  library(dbarts)
  
  X <- data[, -which(names(data) %in% c("Z", "Y"))]  # Covariates
  Z <- data$Z  # Treatment variable
  Y <- data$Y  # Outcome variable
  aerf <- NULL
  
  n <- nrow(X)
  X_df <- as.data.frame(X)
  
  # Generalized Propensity Score (GPS) estimation
  gps_model <- SuperLearner(Y = Z, X = X_df, family = gaussian(), SL.library = "SL.xgboost")
  mu_zx <- predict(gps_model, newdata = X_df)$pred
  gps <- dnorm(Z, mean = mu_zx, sd = sd(Z - mu_zx))

  # Step 1: Estimate IERF using chosen method
  if (method == "superlearner") {
    outcome_model <- SuperLearner(Y = Y, X = data.frame(Z = Z, X_df), family = gaussian(), SL.library = "SL.xgboost")
    mu_hat_z_i <- predict(outcome_model, newdata = data.frame(Z = Z, X_df))$pred
    mu_avg_z_i <- sapply(1:n, function(i) {
      mean(predict(outcome_model, newdata = data.frame(Z = Z[i], X_df))$pred)
    })

    if (!is.null(z_grid)) {
      aerf <- sapply(z_grid, function(z) {
        mean(predict(outcome_model, newdata = data.frame(Z = rep(z, n), X_df))$pred)
      })
    }

  } else if (method == "spline") {
    formula_str <- paste("Y ~ ns(Z, df =", J, ") * (", paste(colnames(X_df), collapse = " + "), ")")
    outcome_model <- lm(as.formula(formula_str), data = data.frame(Y = Y, Z = Z, X_df))
    mu_hat_z_i <- predict(outcome_model)
    mu_avg_z_i <- sapply(1:n, function(i) {
      mean(predict(outcome_model, newdata = data.frame(Z = Z[i], X_df)))
    })

    if (!is.null(z_grid)) {
      aerf <- sapply(z_grid, function(z) {
        mean(predict(outcome_model, newdata = data.frame(Z = rep(z, n), X_df)))
      })
    }

  } else if (method == "doubly_robust") {
    outcome_model <- SuperLearner(Y = Y, X = data.frame(Z = Z, X_df), family = gaussian(), SL.library = "SL.xgboost")
    m_zx <- predict(outcome_model, newdata = data.frame(Z = Z, X_df))$pred
    K_h <- dnorm(0, mean = 0, sd = h) / h  / sum(dnorm(Z - 0, mean = 0, sd = h))
    mu_hat_z_i <- m_zx + (K_h / gps) * (Y - m_zx)
    
    mu_avg_z_i <- sapply(1:n, function(i) {
      m_z_all <- predict(outcome_model, newdata = data.frame(Z = rep(Z[i], n), X_df))$pred
      K_h <- dnorm(Z - Z[i], mean = 0, sd = h) / h / sum(dnorm(Z - Z[i], mean = 0, sd = h)) # Kernel weights for Z[i]
      aipw_vec <- m_z_all + (K_h / gps) * (Y - m_zx)
      mean(aipw_vec)
    })

    if (!is.null(z_grid)) {
      aerf <- sapply(z_grid, function(z) {
        # Compute the average of the doubly robust estimator at z
        m_z_all <- predict(outcome_model, newdata = data.frame(Z = rep(z, n), X_df))$pred
        K_h <- dnorm(Z - z, mean = 0, sd = h) / h / sum(dnorm(Z - z, mean = 0, sd = h)) # Kernel weights for z
        aipw_vec <- m_z_all + (K_h / gps) * (Y - m_zx)
        mean(aipw_vec)
      })
    }

  } else if (method == "dbart") {
    bart_fit <- dbarts::bart(
      x.train = data.frame(Z = Z, X_df),
      y.train = Y,
      keeptrees = TRUE
    )
    mu_hat_z_i <- rowMeans(predict(bart_fit))

    mu_avg_z_i <- sapply(1:n, function(i) {
      X_new <- data.frame(Z = rep(Z[i], n), X_df)
      mean(rowMeans(predict(bart_fit, newdata = X_new)))
    })

    if (!is.null(z_grid)) {
      # Compute the average of the dbart predictions at each z in z_grid
      aerf <- sapply(z_grid, function(z) {
        X_new <- data.frame(Z = rep(z, n), X_df)
        mean(rowMeans(predict(bart_fit, newdata = X_new)))
      })
    }

  } else if (method == "tsbart") {
    fit_tsb <- tsbart::tsbart(
      y = Y,
      tgt = Z,
      x = X_df,
      tpred = Z,
      xpred = X_df,
      nburn = 100,
      nsim = 200,
      ntree = 50,
      use_fscale = TRUE,
      verbose = FALSE
    )
    mu_hat_z_i <- rowMeans(fit_tsb$mcmcdraws)

    mu_avg_z_i <- sapply(1:n, function(i) {
      preds <- predict(fit_tsb, tpred = rep(Z[i], n), xpred = X_df)
      mean(rowMeans(preds))
    })

    if (!is.null(z_grid)) {
      # Compute the average of the tsBART predictions at each z in z_grid
      aerf <- sapply(z_grid, function(z) {
        preds <- predict(fit_tsb, tpred = rep(z, n), xpred = X_df)
        mean(rowMeans(preds))
      })
    }

  } else {
    stop("Invalid method")
  }
  
  diff <- mu_hat_z_i - mu_avg_z_i
  
  # Step 2: Fit tsBART to diff
  fit <- tsbart::tsbart(
    y = diff,
    tgt = Z,
    x = X_df,
    tpred = Z,
    xpred = X_df,
    nburn = nburn,
    nsim = nsim,
    ntree = ntree,
    use_fscale = TRUE,
    verbose = FALSE
  )
  
    # Step 3: Extract rules and compute function norms
  all_rules <- list()
  X_names <- colnames(X_df)
  for (i in seq_along(fit$trees)) {
    for (j in seq_along(fit$trees[[i]])) {
      tree <- fit$trees[[i]][[j]]
      rules <- generate_rule_paths(tree, X_names)
      normalized_rules <- lapply(rules, function(rule) {
        terms <- strsplit(rule, " & ")[[1]]
        paste(sort(terms), collapse = " & ")
      })
      
      for (k in seq_along(normalized_rules)) {
        all_rules <- c(all_rules, normalized_rules[[k]])
      }
    }
  }
  
  all_rules_vec <- unlist(all_rules)
  rule_table <- sort(table(all_rules_vec), decreasing = TRUE)
  rule_set <- names(rule_table)
  total_rules_generated <- length(rule_set)
  rule_table <- as.data.frame(rule_table)
  colnames(rule_table) <- c("rule", "frequency")
  
  # Step 4: Evaluate rule matrix R [n x M]
  R <- sapply(rule_set, function(rule) as.numeric(eval(parse(text = paste0("with(X_df, ", rule, ")")))) )
  R <- as.matrix(R)
  
  # Step 5: Construct spline basis
  Phi <- ns(Z, df = J)  # [n x J]
  
  # Step 6: Build interaction design [n x (M * J)]
  design_list <- lapply(1:ncol(R), function(m) sweep(Phi, 1, R[, m], `*`))
  X_design <- do.call(cbind, design_list)
  groups <- rep(1:ncol(R), each = J)
  
  # Step 7: Estimate GPS (for weights)
  omega <- 1 / gps
  
  # Step 8: Fit group LASSO
  fit_gl <- grpreg(X_design, diff, group = groups, penalty = "grLasso", lambda = lambda, weights = omega)
  selected_groups <- unique(groups[which(coef(fit_gl)[-1] != 0)])
  
  return(list(
    selected_rule_indices = selected_groups,
    selected_rules = rule_set[selected_groups],
    num_selected_rules = length(selected_groups),
    total_rules_generated = total_rules_generated,
    rule_frequencies = rule_table,
    coef = coef(fit_gl),
    gps = gps,  # Original GPS for diagnostics
    diff = diff,
    z_grid = z_grid, 
    aerf = aerf  # Average IERF estimates
  ))
}

generate_rule_paths <- function(tree, X_names) {
  is_leaf <- tree$split_var < 0
  n_nodes <- length(tree$nid)
  
  # Build a lookup of parent → child relationships
  build_tree_edges <- function() {
    parent_ids <- tree$nid[!is_leaf]
    edges <- list()
    for (pid in parent_ids) {
      # tsBART uses nested indexing; children must have parent_id < child_id
      candidates <- tree$nid[tree$nid > pid]
      if (length(candidates) < 2) next
      
      # Assume first two larger nids are left/right
      left <- candidates[1]
      right <- candidates[2]
      
      edges[[as.character(pid)]] <- list(left = left, right = right)
    }
    return(edges)
  }
  
  # Recursively walk from root
  walk_tree <- function(node_id, path) {
    node_idx <- which(tree$nid == node_id)
    if (length(node_idx) == 0) return()
    
    if (is_leaf[node_idx]) {
      rule_paths[[length(rule_paths) + 1]] <<- paste(path, collapse = " & ")
      return()
    }
    
    var_name <- X_names[tree$split_var[node_idx] + 1]
    cut <- round(tree$split_cut[node_idx], 4)
    
    children <- tree_edges[[as.character(node_id)]]
    if (!is.null(children$left)) {
      walk_tree(children$left, c(path, paste0(var_name, " < ", cut)))
    }
    if (!is.null(children$right)) {
      walk_tree(children$right, c(path, paste0(var_name, " >= ", cut)))
    }
  }
  
  rule_paths <- list()
  tree_edges <- build_tree_edges()
  walk_tree(tree$nid[1], character())  # Start from root node
  return(rule_paths)
}

# ---- Inference Step ---- #
inference <- function(data, rule_out, J = 10, z_grid = NULL, x_pred = NULL) {
  library(splines)

  X <- data[, -which(names(data) %in% c("Z", "Y"))]  # Covariates
  Z <- data$Z  # Treatment variable
  Y <- data$Y  # Outcome variable
  n <- nrow(X) 

  if (is.null(x_pred)) x_pred <- X
  n_new <- nrow(x_pred)
  
  diff <- rule_out$diff
  aerf <- rule_out$aerf # Average IERF estimates
  gps <- rule_out$gps  # Generalized propensity scores (for weights)
  rule_set <- rule_out$selected_rules

  if (length(rule_set) == 0) {
    mu_pred_mat <- 
    return(list(
    z_grid = z_grid,
    aerf = aerf,           # average ERF
    rule_set = rule_set,  # selected rules
    rule_matrix = rule_matrix,
    coeffs = C_hat,
    alpha_matrix = alpha_mat,  # rule-wise additive effects αₘ(z)
    mu_pred_mat = mu_pred_mat
  ))
  }

  # Rule matrix R [n x M]
  R_mat <- sapply(rule_set, function(rule) {
    as.numeric(eval(parse(text = rule), envir = as.data.frame(X)))
  })

  M <- ncol(R_mat)
  colnames(R_mat) <- rule_set

  # Spline basis [n x J]
  if (is.null(z_grid)) z_grid <- sort(unique(Z))
  Phi <- ns(Z, df = J)

  # Design matrix X = R ⊗ Phi [n x M*J]
  design_list <- lapply(1:M, function(m) sweep(Phi, 1, R_mat[, m], `*`))
  X_design <- do.call(cbind, design_list)
  groups <- rep(1:M, each = J)

  # Weight matrix
  W <- diag(1/gps)

  # Least squares solution: solve (X^T W X)^{-1} X^T W y
  XtW <- t(X_design) %*% W
  XtWX <- XtW %*% X_design
  C_hat <- solve(XtWX + diag(1e-6, ncol(XtWX)), XtW %*% diff)
  
  # Build alpha_m(z) = sum_j C_{jm} phi_j(z)
  Phi_z <- ns(z_grid, df = J)
  alpha_mat <- matrix(NA, nrow = length(z_grid), ncol = M)
  for (m in 1:M) {
    coeffs <- C_hat[((m - 1) * J + 1):(m * J)]
    alpha_mat[, m] <- Phi_z %*% coeffs
  }
  colnames(alpha_mat) <- rule_set
  
  rule_matrix <- sapply(rule_set, function(rule) {
    as.numeric(eval(parse(text = rule), envir = as.data.frame(x_pred)))
  })
  
  mu_pred_mat <- matrix(aerf, nrow = length(z_grid), ncol = n_new)
  for (j in 1:n_new) {
    mu_pred_mat[, j] <- aerf + alpha_mat %*% rule_matrix[j, ]
  }
  
  return(list(
    z_grid = z_grid,
    aerf = aerf,           # average ERF
    rule_set = rule_set,  # selected rules
    rule_matrix = rule_matrix,
    coeffs = C_hat,
    alpha_matrix = alpha_mat,  # rule-wise additive effects αₘ(z)
    mu_pred_mat = mu_pred_mat
  ))
}

estimate_marginal_density <- function(z_grid, n = 1000, phi = 10) {
  X <- simulate_X(n)
  pi_vals <- pi_x(X)
  density_matrix <- sapply(z_grid, function(z) {
    dbeta(z, shape1 = phi * pi_vals, shape2 = phi * (1 - pi_vals))
  })
  density_means <- colMeans(density_matrix)
  density_means / sum(density_means)
}

# calc_bias_rmse <- function(S = 50, 
#                            n = 500,
#                            z_grid = seq(0.1, 0.9, length.out = 100),
#                            phi = 10,
#                            sigma = 1,
#                            k = 1,
#                            estimation_fn) {
#   L <- length(z_grid)
  
#   # Precompute marginal density weights ω̄(z)
#   w_z <- estimate_marginal_density(z_grid = z_grid, n = 10000, phi = phi)
  
#   # Fix X across simulations
#   X <- simulate_X(n)
  
#   # Compute true μ_i(z) matrix [n x L]
#   f_x <- X[, 1]  + X[, 3] + X[, 4]
#   r1 <- as.numeric(X[, 1] == 1 & X[, 2] == 0)
#   r2 <- as.numeric(X[, 5] == 1 & X[, 6] == 0)
  
#   mu_true_matrix <- outer(
#     z_grid,
#     1:n,
#     Vectorize(function(z, i) f_x[i] + k * (z * r1[i] - z * r2[i]))
#   )
#   mu_true_matrix <- t(mu_true_matrix)  # [n x L]
  
#   # Store simulation results: [n x L x S]
#   mu_hat_array <- array(NA, dim = c(n, L, S))
  
#   for (s in 1:S) {
#     sim_data <- generate_sim_data(X = X, phi = phi, sigma = sigma, k = k)
#     Z <- sim_data$Z
#     Y <- sim_data$Y
    
#     mu_hat_matrix <- estimation_fn(z_grid = z_grid, X = X, Z = Z, Y = Y)
#     mu_hat_array[, , s] <- mu_hat_matrix
#     print(s)
#   }
  
#   # Compute bias and RMSE vectors per individual
#   bias_vec <- numeric(n)
#   rmse_vec <- numeric(n)
  
#   for (i in 1:n) {
#     mu_true <- mu_true_matrix[i, ]
#     mu_hat_mean <- rowMeans(mu_hat_array[i, , ])  # avg across S
    
#     # Bias: ∑_z |mean_hat - true| * w(z)
#     bias_vec[i] <- sum(abs(mu_hat_mean - mu_true) * w_z)
    
#     # RMSE: sqrt of ∑_z [mean (hat - true)^2] * w(z)
#     mse_z <- sapply(1:L, function(j) {
#       mean((mu_hat_array[i, j, ] - mu_true_matrix[i, j])^2)
#     })
#     rmse_vec[i] <- sqrt(sum(mse_z * w_z))
#   }
  
#   data.frame(bias = bias_vec, rmse = rmse_vec)
# }


########### Simulation
# Rule Discovery
simulate_rule_discovery <- function(methods = c("superlearner"),
                                    k_vals = 1:10,
                                    S = 10,
                                    n = 500,
                                    lambda = 0.1,
                                    h = 0.05,
                                    nburn = 50,
                                    nsim = 10,
                                    ntree = 40,
                                    J = 10,
                                    phi = 1, 
                                    sigma = 0.5,
                                    true_rules = c("X5 >= 1 & X6 < 1", "X1 >= 1 & X2 < 1")) {
  results <- list()

  for (method in methods) {
    for (k in k_vals) {
      cat("Running method:", method, "k =", k, "\n")
      tp_list <- fp_list <- fn_list <- numeric(S)

      for (s in 1:S) {
        X <- simulate_X(n = n)
        data <- generate_sim_data(X = X, phi = phi, sigma = sigma, k = k)
        out <- rule_generation(data, J = J, lambda = lambda, h = h, nburn = nburn, nsim = nsim, ntree = ntree, method = method)

        if (length(out$selected_rules) == 0) {
          tp <- 0; fp <- 0; fn <- length(true_rules)
        } else {
          selected <- out$selected_rules
          tp <- sum(selected %in% true_rules)
          fp <- sum(!selected %in% true_rules)
          fn <- sum(!true_rules %in% selected)
        }

        tp_list[s] <- tp
        fp_list[s] <- fp
        fn_list[s] <- fn
      }

      recall <- mean(tp_list / (tp_list + fn_list + 1e-10))
      precision <- mean(tp_list / (tp_list + fp_list + 1e-10))
      f1_score <- mean(2 * (precision * recall) / (precision + recall + 1e-10))

      results[[length(results) + 1]] <- data.frame(
        method = method,
        k = k,
        recall = recall,
        precision = precision,
        f1_score = f1_score
      )
    }
  }

  return(do.call(rbind, results))
}

simulate_inference <- function(methods = c("superlearner", "spline", "aipw", "dbart", "tsbart"),
                                    k = 3, n = 100, S = 10, J = 10, lambda = 0.1, h = 0.05, nburn = 50, nsim = 10, ntree = 40,
                                    z_grid = seq(0.1, 0.9, length.out = 100),
                                    phi = 1, sigma = 0.1) {
    L <- length(z_grid)
    w_z <- estimate_marginal_density(z_grid = z_grid, n = 10000, phi = phi)
    X <- simulate_X(n)

    f_x <- X[, 1]  + X[, 3] + X[, 4]
    f_z <- 1 + Z + Z^2 - Z^3
    r1 <- as.numeric(X[, 1] == 1 & X[, 2] == 0)
    r2 <- as.numeric(X[, 5] == 1 & X[, 6] == 0)

    mu_true_matrix <- t(outer(
    z_grid,
    1:n,
    Vectorize(function(z, i) f_z[i] + 2 * k * (z * r1[i] - z * r2[i]))
    ))  

    results <- list()

    for (method in methods) {
        cat("Method:", method, "\n")
        cre_bias_vec <- numeric(S)
        cre_rmse_vec <- numeric(S)
        bias_vec <- numeric(S)
        rmse_vec <- numeric(S)

        for (s in 1:S) {
            sim_data <- generate_sim_data(X = X, phi = phi, sigma = sigma, k = k)

            # CRE
            rule_out <- rule_generation(sim_data, z_grid = z_grid, J = J, h = h, nburn = nburn, nsim = nsim, ntree = ntree, lambda = lambda, method = method)
            gps <- rule_out$gps  # Get the gps from rule generation
            inf <- inference(sim_data, rule_out, J = J, z_grid = z_grid, x_pred = NULL)
            mu_hat_matrix <- t(inf$mu_pred_mat)
            mu_diff <- mu_hat_matrix - mu_true_matrix
            cre_bias_vec[s] <- sum(abs(colMeans(mu_diff)) * w_z)
            cre_rmse_vec[s] <- sum(sqrt(colMeans(mu_diff^2)) * w_z)

            # Baseline
            X <- sim_data[, -which(names(data) %in% c("Z", "Y"))]  # Covariates
            Z <- sim_data$Z  # Treatment variable
            Y <- sim_data$Y  # Outcome variable
            X_df <- as.data.frame(X)
            mu_hat_matrix <- matrix(NA, nrow = n, ncol = L)

            if (method == "superlearner"){
                outcome_model <- SuperLearner(Y = Y, X = data.frame(Z = Z, X_df), family = gaussian(), SL.library = "SL.xgboost")
                for (j in 1:L) {
                    z_val <- z_grid[j]
                    new_data <- data.frame(Z = rep(z_val, n), X_df)
                    mu_hat_matrix[, j] <- predict(outcome_model, newdata = new_data)$pred
                }
            } else if (method == "spline"){
                formula_str <- paste("Y ~ ns(Z, df =", J, ") * (", paste(colnames(X_df), collapse = " + "), ")")
                outcome_model <- lm(as.formula(formula_str), data = data.frame(Y = Y, Z = Z, X_df))
                for (j in 1:L) {
                    z_val <- z_grid[j]
                    new_data <- data.frame(Z = rep(z_val, n), X_df)
                    mu_hat_matrix[, j] <- predict(outcome_model, newdata = new_data)
                }
            } else if (method == "doubly_robust") {
                # Fit a model to estimate the outcome
                outcome_model <- SuperLearner(Y = Y, X = data.frame(Z = Z, X_df), family = gaussian(), SL.library = "SL.xgboost")
                m_zx <- predict(outcome_model, newdata = data.frame(Z = Z, X_df))$pred
                for (j in 1:L) {
                    z_val <- z_grid[j]
                    new_data <- data.frame(Z = rep(z_val, n), X_df)
                    m_z_all <- predict(outcome_model, newdata = new_data)$pred
                    K_h <- dnorm(Z - z_val, mean = 0, sd = h) / h / sum(dnorm(Z - z_val, mean = 0, sd = h))
                    mu_hat_matrix[, j] <- m_z_all + (K_h / gps) * (Y - m_zx)
                }
            } else if (method == "dbart") {
                outcome_model <- dbarts::bart(x.train = data.frame(Z = Z, X_df), y.train = Y, keeptrees = TRUE)
                for (j in 1:L) {
                    z_val <- z_grid[j]
                    new_data <- data.frame(Z = rep(z_val, n), X_df)
                    mu_hat_matrix[, j] <- rowMeans(predict(outcome_model, newdata = new_data))
                }
            } else if (method == "tsbart") {
                outcome_model <- tsbart::tsbart(
                    y = Y,
                    tgt = Z,
                    x = X_df,
                    tpred = Z,
                    xpred = X_df,
                    nburn = 100,
                    nsim = 200,
                    ntree = 50,
                    use_fscale = TRUE,
                    verbose = FALSE
                    )
                for (j in 1:L) {
                    z_val <- z_grid[j]
                    preds <- predict(outcome_model, tpred = rep(z_val, n), xpred = X_df)
                    mu_hat_matrix[, j] <- rowMeans(preds)
                }
            }

            mu_diff <- mu_hat_matrix - mu_true_matrix
            bias_vec[s] <- sum(abs(colMeans(mu_diff)) * w_z)
            rmse_vec[s] <- sum(sqrt(colMeans(mu_diff^2)) * w_z)
        }

        results[[paste0("CRE_", method)]] <- data.frame(
        bias_mean = mean(cre_bias_vec),
        bias_sd = sd(cre_bias_vec),
        rmse_mean = mean(cre_rmse_vec),
        rmse_sd = sd(cre_rmse_vec)
        )

        results[[method]] <- data.frame(
        bias_mean = mean(bias_vec),
        bias_sd = sd(bias_vec),
        rmse_mean = mean(rmse_vec),
        rmse_sd = sd(rmse_vec)
        )
    }

    do.call(rbind, results)
}

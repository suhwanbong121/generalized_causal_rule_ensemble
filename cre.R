X_df <- as.data.frame(X)
X_input <- data.frame(Z = Z, X_df)

sl_fit <- SuperLearner(
  Y = Y,
  X = X_input,
  family = gaussian(),
  SL.library = SL.library
)

X_input_new <- data.frame(Z = Z, X_df)

# Predict at observed (Z_i, X_i)
mu_hat <- predict(sl_fit, newdata = X_input)$pred
mu_hat_bar <- rep(NA, 100)
for (i in 1:100){
  mu_hat_bar <- predict(sl_fit, newdata = data.frame(Z = Z[i], X_df))$pred
}
residuals <- mu_hat - mu_hat_bar

mean(residuals)

# 1. Fit tsBART
fit <- tsbart::tsbart(
  y = residuals,
  tgt = Z,
  x = as.data.frame(X),
  verbose = FALSE,
  nburn = 50,
  nsim = 200,
  ntree = 20,
  base_tree = 0.99,      # more aggressive growth
  power_tree = 1.5,
  sd_control = 5,
  ecross = 0.2,
  use_fscale = TRUE
)

ntree = 20,            # fewer trees → each tree learns more
base_tree = 0.99,      # more aggressive growth
power_tree = 1.5
plot(Z[order(Z)], mu_vec[order(Z)], type = "l", lwd = 2)

Z[order(Z)]

diff(sort(unique(Z)))
mu_vec <- trees[[1]][[1]]$mu[[3]]
tree1 <- trees[[1]][[1]]
data.frame(
  node = tree1$nid,
  split_var = tree1$split_var,
  split_cut = tree1$split_cut,
  mu_length = lengths(tree1$mu)
)

plot(sort(unique(Z)), trees[[1]][[1]]$mu[[3]], type = "l", col = "red")
lines(sort(unique(Z)), trees[[2]][[1]]$mu[[3]], col = "blue")
lines(sort(unique(Z)), trees[[3]][[1]]$mu[[3]], col = "green")
lines(sort(unique(Z)), trees[[100]][[1]]$mu[[3]], col = "black")
sapply(1:5, function(i) trees[[i]][[1]]$split_var)

plot_node_mu <- function(fit, iter = 1, tree = 1, node = 1, z_vals = NULL, col = "blue", lwd = 2, ...) {
  mu_vec <- fit$trees[[iter]][[tree]]$mu[[node]]
  
  if (is.null(z_vals)) {
    # Reconstruct tref (assumed to be sort(unique(Z)))
    z_vals <- seq(0, 1, length.out = length(mu_vec))  # or use actual Z grid if known
  }
  
  plot(z_vals, mu_vec, type = "l", col = col, lwd = lwd,
       xlab = "z", ylab = expression(mu(z)), 
       main = paste0("Tree ", tree, ", Node ", node, ", Iter ", iter),
       ...)
}

par(mfrow = c(1,1))
for (i in 1:4) {
  plot_node_mu(fit, iter = 1, tree = 1, node = i)
}

extract_rules_from_tree <- function(tree, X_names) {
  n_nodes <- length(tree$nid)
  rules <- vector("list", n_nodes)
  
  for (i in seq_len(n_nodes)) {
    if (tree$split_var[i] < 0) next  # leaf node
    
    var_idx <- tree$split_var[i] + 1  # 0-based to 1-based
    var_name <- X_names[var_idx]
    cut <- tree$split_cut[i]
    
    rule_str <- paste0(var_name, " <= ", round(cut, 3))
    rules[[i]] <- rule_str
  }
  
  rules <- rules[!sapply(rules, is.null)]
  return(unique(unlist(rules)))
}

tree1 <- fit$trees[[300]][[43]]
X_names <- colnames(X)
rules <- extract_rules_from_tree(tree1, X_names)
print(rules)

rule_matrix <- sapply(rules, function(rule) {
  as.numeric(eval(parse(text = paste0("with(X, ", rule, ")"))))
})

all_rules <- list()
for (i in 1:nsim) {
  for (j in 1:ntree) {
    rules <- extract_rules_from_tree(fit$trees[[i]][[j]], X_names)
    all_rules <- c(all_rules, rules)
  }
}
unique_rules <- unique(unlist(all_rules))

tree <- fit$trees[[100]][[1]]
internal_nodes <- which(tree$split_var >= 0)
terminal_nodes <- setdiff(seq_along(tree$nid), internal_nodes)
mu_vec <- tree$mu[[terminal_nodes[1]]]  # This is μ_j(z) for one leaf


rules <- extract_rules(treeList, X, max_depth = 3)
plot(sort(unique(Z)), trees[[1]][[1]]$mu[[3]])
trees[[1]][[1]]$mu[[1]]
extract_rules_from_tsbart <- function(fit, X) {
  n <- nrow(X)
  trees <- fit$trees
  M <- length(trees)
  rule_matrix <- matrix(0, nrow = n, ncol = M)
  
  for (m in seq_len(M)) {
    tree <- trees[[m]]
    
    # predict with only tree m
    pred_m <- predict(fit, xpred = X, treeind = m, alltrees = FALSE)
    
    # Create a binary rule indicator using a threshold (e.g., median)
    rule_matrix[, m] <- as.numeric(pred_m > median(pred_m))
  }
  
  colnames(rule_matrix) <- paste0("r", seq_len(M))
  return(rule_matrix)
}


# 2. Extract rule matrix R_mat [n x M]
R_mat <- extract_rules_from_tsbart(fit, X)

# 3. Spline basis for Z
library(splines)
B_z <- ns(Z, df = 5)

# 4. Construct interaction design matrix [n x (M * J)]
design_list <- lapply(1:ncol(R_mat), function(m) {
  sweep(B_z, 1, R_mat[, m], `*`)
})
X_group_lasso <- do.call(cbind, design_list)

# 5. Fit group lasso
library(grpreg)
groups <- rep(1:ncol(R_mat), each = ncol(B_z))
fit_gl <- grpreg(X_group_lasso, residuals, group = groups, penalty = "grLasso")

# 6. Extract selected rule indices
selected_rules <- unique(groups[which(coef(fit_gl)[-1] != 0)])







tree$split_var
tree <- fit$trees[[50]][[17]]
X_names <- colnames(X)

rules <- generate_rule_paths(tree, X_names)
print(rules)


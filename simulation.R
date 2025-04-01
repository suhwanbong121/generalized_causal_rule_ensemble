library(SuperLearner)
library(splines)
library(grpreg)
library(tsbart)
library(dbarts)
library(ggplot2)

source("functions.R")

set.seed(0)
# methods = c("superlearner", "spline", "doubly_robust", "dbart", "tsbart"),
# 1. Rule generation
true_rules <- c("X5 >= 1 & X6 < 1", "X1 >= 1 & X2 < 1")
sim <- simulate_rule_discovery(
    methods = c("superlearner", "spline", "doubly_robust"),
    k_vals = 1:5,
    S = 20,
    n = 200,
    lambda = 0.5,
    h = 0.05,
    nburn = 50,
    nsim = 10,
    ntree = 40,
    J = 5,
    phi = 1, 
    sigma = 0.1,
    true_rules = true_rules
    )

sim

## Effect modifier ..?

# Save Recall plot
recall_plot <- ggplot(sim, aes(x = k, y = recall, color = method)) +
  geom_line() + geom_point() +
  ylab("Recall") + ggtitle("Recall by method and k")
ggsave("figures/recall_plot.png", plot = recall_plot, width = 8, height = 6)

# Save Precision plot
precision_plot <- ggplot(sim, aes(x = k, y = precision, color = method)) +
  geom_line() + geom_point() +
  ylab("Precision") + ggtitle("Precision by method and k")
ggsave("figures/precision_plot.png", plot = precision_plot, width = 8, height = 6)

# Save F1 Score plot
f1_plot <- ggplot(sim, aes(x = k, y = f1_score, color = method)) +
  geom_line() + geom_point() +
  ylab("F1 score") + ggtitle("F1 score by method and k")
ggsave("figures/f1_score_plot.png", plot = f1_plot, width = 8, height = 6)

# 2. Estimation
res <- simulate_inference(methods = c("superlearner", "spline", "doubly_robust"),
    k = 5,
    n = 300,
    S = 5,
    J = 5,
    phi = 1,
    sigma = 0.1,
    nburn = 50,
    nsim  = 10,
    ntree = 50,
    z_grid = seq(0.1, 0.9, length.out = 100),
    lambda = 0.3)

res

write.csv(res, "data_output.csv", row.names = FALSE)


n = 200
lambda = 1
k = 5
X = simulate_X(n)
data = generate_sim_data(X, phi, sigma, k = k)
rule <- rule_generation(data, J = 20, lambda = 0.3, method = "superlearner")

rule$selected_rules

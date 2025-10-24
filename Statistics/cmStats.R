import("modules")

# cmStats module definition: 
cmStats <- module({
import("caracas")
import("stats")
  
  .foo <- identity
  caracas::def_sym("Evn", "x", "X", "y", "Y", "z", "n", "i", "S", "N", "p", "k", "y_hat", "Y_bar", "X_bar", "X_hat", "n_goup", "X_goup", "Lambda", "e", "sigma", "beta", "theta", "Theta", "Theta_hat", "sigma_x", "sigma2_x", "Sigma", "mu", "mu_x", "mu_0", "rho", "sum_x", "phi", "alpha", "Z_alpha")
 
  # sympy functions for symbolic computation:
  
  factorial_sym  <- \(x) caracas::sympy_func(x, "factorial")
  sum_sym <- caracas::sum_
  sqrt_sym <-  \(x) caracas::sympy_func(x, "sqrt")
  min_sym <- \(x) caracas::sympy_func(x, "Min")
  max_sym <-\(x) caracas::sympy_func(x, "Max")

  # major statistics symbolic functions:
  mu_sym <- sum_sym(X, X, i, n, doit = FALSE) / N
  P_sym <- Evn/p
  Z_score_sym <- \(X=X) (X - mu) / sigma
  Z_sym <- (X - mu) / (sigma / sqrt_sym(n))
  Z_small_sym <- (X - mu) / (S / sqrt_sym(n))
  
  rho_sym <- \(X=X, Y=Y) sum(Z_score_sym(X) * Z_score_sym(Y))/ (n - 1)
  σ_2_sym <- (sum_x - mu)^2/N
  variance_sym <- sum(X - mu_sym)^2/ (n - x)
  SS_B_sym <-  sum_(n_goup * (X_goup - X)^2 , i, 1, n, doit = FALSE)
  SS_sym <- sum_((x - X )^2, i, 1, n, doit = FALSE)
  S_2_sym <- sum_((x - X )^2 / (n - 1), i, 1, n, doit = FALSE)
  SS_W_sym <- SS_sym
  SS_T_sym <- SS_sym
  C_n_k_sym <- factorial(n) / (factorial(k) * factorial(n - k))
  P_n_k_sym <- n/(n - k)
  μ_hat_sym <- min_sym(X) + max_sym(X)/2
  odds_sym <- p/(1 - p)
  logit_sym <-  log(p/(1- p))
  
  P <-  function(Evn, p) caracas::as_func(P_sym)(Evn, p)
  
  bayes_rule_sym <- (P(X, Theta) * P(Theta, y)) / P(X, y)
  abs_loss_sym <- abs(Theta - Theta_hat)
  sqrt_loss_sym <- (Theta - Theta_hat)^2
  
  
  # test new functions implementation as symbolic
  #caracas::sum_(X, X, i, n)/N
  #caracas::sympy_func(y, "factorial")
  #factorial_sym(n) / (factorial_sym(k) * factorial_sym((n - k)))
  #as_func(mu_sym)

  
  .Xs <-  c(1:10)
  
  Σ <-  function(xs, n=length(xs)) sum(xs, n)
  Σ_sym <- as_func(caracas::sum_(x, x, i, n))

  μ <-  function(xs, N = length(xs)) {
    sym_fn <- caracas::as_func(mu_sym); sym_fn(N, Σ(xs, 0))
  }
  
  σ_2 <-  function (xs) Σ((xs - μ(xs))^2, length(xs)) / length(xs)
  
  σ_2_test <-  function (xs) Σ((xs - μ(xs))^2, length(xs)) / length(xs)
  
  
  σ_2_sample <-  function (xs) Σ((xs - μ(xs))^2, length(xs))/ (length(xs) - 1)
  
  ρ <- function(Xs, Ys) {
    caracas::as_func(rho_sym(X, Y))(stats::mean(Xs), length(Xs), sigma=stats::sd(Xs), Xs, Ys)
    
  }
  
  f <-  function(x=x) x
  
  Z_score <- function(Xs=Xs) {
    caracas::as_func(Z_score_sym(X))(mu = mean(Xs), sigma = stats::sd(Xs), Xs)
  }
  
  Z_fn <- function(mu_0 = 0, n=length(Xs), sigma = stats::sd(Xs), Xs = mean(Xs)){
    caracas::as_func(cmStats$Z_sym)(mu_0, n, sigma, Xs)
  }
  
  C_n_k <- function(n, k) {
    caracas::as_func(C_n_k_sym)(k, n)
  }
  
  P_n_k <- function(n, k) {
    caracas::as_func(P_n_k_sym)(k, n)
  }
  
})


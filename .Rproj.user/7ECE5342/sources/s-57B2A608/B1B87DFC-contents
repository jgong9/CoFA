#' A function to estimate the variance of random noise in multivariate data
#' @export var_omega_est
#' @param z_data_cent The matrix of centered multivariate data
#' @param V0_mat The matrix of estimated eigen vectors for common components
#' @param est_beta_vec The vector of the estimated beta values
#' @param est_eta_vec The vector of the estimated eta values
#' @param est_K0 The estimated K0 value
#' @return An object of a solution from CVXR R package
#' @import CVXR

var_omega_est <- function(z_data_cent, V0_mat, est_beta_vec, est_eta_vec, est_K0){

  n <- dim(z_data_cent)[1]
  Cov_Z <- t(z_data_cent) %*% z_data_cent / (n-1)
  P_v0 <- V0_mat  %*% solve( t(V0_mat) %*% V0_mat  )%*% t(V0_mat)

  ## estimation for sigma^{2}_{omega}
  ## based on Ordinary Least Square with a constraint
  if(est_K0 == 1){
    product_hat <- P_v0 %*% Cov_Z + Cov_Z %*% P_v0 - P_v0 %*% Cov_Z %*% P_v0 -  est_beta_vec[1] * est_eta_vec[1] *  V0_mat  %*% t(V0_mat)
  } else{
    product_hat <- P_v0 %*% Cov_Z + Cov_Z %*% P_v0 - P_v0 %*% Cov_Z %*% P_v0 -
    V0_mat %*% diag(  est_beta_vec[1:est_K0] * est_eta_vec[1:est_K0]) %*% t(V0_mat)
  }
  y_vec <- as.numeric(product_hat)
  design_mat <- as.numeric(P_v0)
  betaHat <- Variable(1)
  objective <- Minimize(sum((y_vec - design_mat %*% betaHat)^2))

  problem <- Problem(objective, constraints = list(betaHat >= 0))

  solution_ols <- solve(problem)
  res <- solution_ols$getValue(betaHat)
  return(res)
  }

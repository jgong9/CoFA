#' Common Factor Analysis between Functional and Multivariate Data by Adaptive Nuclear Norm Penalty
#'
#' @export CoFA
#' @param y_sample The functional data matrix
#' @param z_sample The multivariate data matrix
#' @param y_center A logical value to indicate whether y_sample needs to be be centered (TRUE) or not (FALSE)
#' @param z_standard A logical value to indicate whether z_sample needs to be standardized (TRUE) or not (FALSE)
#' @param t_vec The observation points of functional data
#' @param num_knot An integer to choose the total numbers of knots
#' @param user_K1 An integer for K1 fixed by users. If not specified (NULL), estimation by BEMA method applies
#' @param user_K2 An integer for K2 fixed by users. If not specified (NULL), estimation by BEMA method applies
#' @param gamma The parameter used for weights in adaptive nuclear norm penalty
#' @param BEMA_alpha One of the two tuning parameters in the BEMA method. It determines a percentage of nonzero eigenvalue in the middle range. Default is 0.2
#' @param BEMA_beta One of the two tuning parameters in the BEMA method. It controls the probability of over-estimating K1 and K2. Recommended default is 0.1
#' @param control A list that contains the cross validation setting
#' @return A list that contains all estimation and prediction results
#' @importFrom face select.knots pspline
#' @importFrom splines splineDesign
#' @importFrom RMTstat qmp ptw

CoFA <- function(y_sample, z_sample, t_vec, y_center=TRUE, z_standard=TRUE, num_knot,
                 num_K1=NULL, num_K2=NULL, gamma=2, BEMA_alpha = 0.2, BEMA_beta = 0.1,
                 control = list(num_rep=100, num_fold=5, num_tau=100)
){

  ### 1. Data pre-processing
  n <- dim(y_sample)[1]
  p <- dim(z_sample)[2]
  m <- length(t_vec)

  ## 1-a. Multivariate data standardization
  if(z_standard==TRUE){
    z_sample_std <- scale(z_sample)
  } else {
    z_sample_std <- z_sample
  }

  ## 1-b. Functional data centering
  # Based on cubic P-spline smoothing
  if(num_knot < 12){
    stop("Please choose the number of knots at least greater than 12, which will have 4+4 outer knots and 4 inner knots.")
  }
  knots_set <- (num_knot - 7)
  knots_vec <- select.knots(t_vec, knots = knots_set ,p=3)
  B_mat <- splineDesign(knots_vec, t_vec, ord = 4)

  w <- quadWeights(argvals = t_vec , method = "trapezoidal")
  W_mat <- diag(w)
  W_sqrt <- diag( sqrt(w) )
  G <- t(B_mat) %*% W_mat %*% B_mat
  G_sqrt <- eigen(G)$vectors %*% diag( sqrt(eigen(G)$values)  ) %*% t( eigen(G)$vectors )
  G_sqrt_inverse <- eigen(G)$vectors %*% diag( 1/sqrt(eigen(G)$values)  ) %*% t( eigen(G)$vectors )
  B_tilde <- B_mat %*% G_sqrt_inverse

  # Demeaning functional data
  if(y_center==TRUE){
    id <- rep(1:n,each=m)
    t_pspline <- rep(t_vec,times=n)
    y_stack <- as.vector(t(y_sample))
    data <- data.frame(y=y_stack,
                       argvals = t_pspline,
                       subj = id)
    fit_pspline <- pspline(data,knots=knots_set)
    y_centered <- t(matrix(fit_pspline$fitted.values, length(t_vec), n))
    y_demeaned <- y_sample - y_centered
    y_demeaned <- as.matrix(y_demeaned)
  } else {
    y_centered <- t(matrix(0, length(t_vec), n))
    y_demeaned <- y_sample
  }


  ### 2. CoFA for common component part
  ## 2-a. Estimation for K0, phi_k0, and v_k0
  CoFA_CrossCov_fit <- CoFA_CrossCov(Y_cent = y_demeaned, Z_cent = z_sample_std, B_tilde = B_tilde, t_vec=t_vec, gamma = gamma,
                            control = list(num_rep= control$num_rep, num_fold=control$num_fold, num_tau=control$num_tau))

  ## 2-b. Smooth covariance function estimation by fpca.sc of refund R package
  fpca_sc_fit <- fpca_sc_modified(Y = y_demeaned, ydata = NULL, Y.pred = NULL, argvals = t_vec, random.int = FALSE,
                                  nbasis = dim(B_tilde)[2] , pve = 0.99, npc = NULL, var = TRUE, simul = FALSE, sim.alpha = 0.95,
                                  useSymm = FALSE, makePD = FALSE, center = FALSE, cov.est.method = 2, integration = "trapezoidal")
  cov_fpcaSC_cov_hat <- fpca_sc_fit$cov.hat
  est_error_epsilon_fpca_sc <- fpca_sc_fit$sigma2

  ## 2-c. Estimation for all parameters related to common factors
  est_K0 <- CoFA_CrossCov_fit$rank

  est_eta <- vector(mode="numeric", length = est_K0)
  Cross_cov <- 1/(n-1) *( t(y_demeaned) %*% (z_sample_std) )
  Gamma_mat <- t(B_tilde) %*% W_mat %*% Cross_cov

  est_lambda0 <- vector(mode="numeric", length = est_K0)
  est_beta <- vector(mode="numeric", length = est_K0)
  for(j in 1:est_K0){
    est_eta[j] <- t(CoFA_CrossCov_fit$U[,j]) %*% Gamma_mat %*% CoFA_CrossCov_fit$V[,j]
    est_beta[j] <- est_eta[j] / ( t(B_tilde %*% CoFA_CrossCov_fit$U[,j]) %*%  W_mat %*% cov_fpcaSC_cov_hat %*%  W_mat %*% B_tilde %*% CoFA_CrossCov_fit$U[,j] )
    est_lambda0[j] <-   (t(B_tilde %*% CoFA_CrossCov_fit$U[,j]) %*%  W_mat %*% cov_fpcaSC_cov_hat %*%  W_mat %*% B_tilde %*% CoFA_CrossCov_fit$U[,j]  )
  }



  ### 3. Estimation for independent component part
  ## 3-a functional data
  # Estimatiion K_{1} by BEMA and MLE of sigma^{2}_{epsilon}
  cov_y <- t(y_demeaned) %*% (y_demeaned) / (n-1)

  svd_Cov_y <- svd(cov_y)
  d_Cov_y <- svd_Cov_y$d

  p_tilde_fnt <- min(n, m)
  gamma_n_fnt <- m/n
  bulk_eigenvalue_order_fnt <- seq( ceiling(p_tilde_fnt * BEMA_alpha), floor( p_tilde_fnt * (1-BEMA_alpha) ) )
  bulk_eigenvalue_fnt <- d_Cov_y[bulk_eigenvalue_order_fnt]
  q_k_vec_fnt <-  qmp( bulk_eigenvalue_order_fnt/p_tilde_fnt , ndf=n, pdim=m, var=1, lower.tail = FALSE, log.p = FALSE )

  est_error_epsilon_BEMA <- sum( bulk_eigenvalue_fnt * q_k_vec_fnt ) / sum(q_k_vec_fnt^(2))

  t_1minusbeta <- ptw( 1 - BEMA_beta , beta=1, lower.tail = TRUE, log.p = FALSE)
  T_hat_fnt <- est_error_epsilon_BEMA *
    ( (1+sqrt(gamma_n_fnt))^(2) + t_1minusbeta * n ^(-2/3) * (gamma_n_fnt)^(-1/6) * (1+sqrt(gamma_n_fnt))^(4/3) )

  est_K0_plus_K1 <- sum( d_Cov_y > T_hat_fnt )
  est_error_epsilon_MLE  <-   mean( d_Cov_y[-c(1:est_K0_plus_K1)] )
  est_K1 <- est_K0_plus_K1 - est_K0
  if(is.null(num_K1) == TRUE){
    if(est_K1 <= 0){
      cat("Estimation results in K0 >= K0+K1. K1 is set to 3.\n ")
      est_K1 <- 3
    }
  } else {
    est_K1 <- num_K1
  }

  # Estimatiion for phi_k1 and lambda_k1
  est_lambda1 <- vector(mode="numeric", length = est_K1)
  U0 <- CoFA_CrossCov_fit$U[,1:est_K0]
  P_U1 <- diag(1, dim(B_tilde)[2]) - U0 %*% solve( t(U0) %*% U0) %*% t(U0)
  svd_indep_fun <- svd( P_U1 %*% t(B_tilde) %*% ( W_mat %*% cov_y %*% W_mat  ) %*% B_tilde  %*% t(P_U1) )
  if(est_K0_plus_K1 > dim(P_U1)[1]){
    stop(
      cat("The estimated number of all factors in functional data is greater than the number of basis functions.\n",
            "Choose \"num_knot\" >=", est_K0_plus_K1+4, " or specify \"user_K1\" argument.\n"
            )
    )
  }



  for(j in 1:min(c(est_K1, length(svd_indep_fun$d) ) ) ){
    est_lambda1[j] <-   (t(B_tilde %*% svd_indep_fun$u[,j]) %*% W_mat %*% cov_fpcaSC_cov_hat %*% W_mat %*% B_tilde %*% svd_indep_fun$u[,j]  )
  }



  ## 3-b multivariate data
  # Estimatiion K_{2} by BEMA and MLE of sigma^{2}_{omega}
  Cov_Z <- 1/(n-1) * t(z_sample_std) %*% z_sample_std
  svd_Cov_z <- svd(Cov_Z)
  d_Cov_z <- svd_Cov_z$d

  p_tilde <- min(n, p)
  gamma_n <- p/n
  bulk_eigenvalue_order <- seq( ceiling(p_tilde * BEMA_alpha), floor( p_tilde * (1-BEMA_alpha) ) )
  bulk_eigenvalue <- d_Cov_z[bulk_eigenvalue_order]
  q_k_vec <-  qmp( bulk_eigenvalue_order/p_tilde , ndf=n, pdim=p, var=1, lower.tail = FALSE, log.p = FALSE )

  est_error_omega_BEMA <- sum( bulk_eigenvalue * q_k_vec ) / sum(q_k_vec^(2))

  t_1minusbeta <- ptw( 1 - BEMA_beta , beta=1, lower.tail = TRUE, log.p = FALSE)
  T_hat <- est_error_omega_BEMA *
    ( (1+sqrt(gamma_n))^(2) + t_1minusbeta * n ^(-2/3) * (gamma_n)^(-1/6) * (1+sqrt(gamma_n))^(4/3) )

  est_K0_plus_K2 <- sum( d_Cov_z > T_hat )
  est_error_omega_MLE  <-   mean( d_Cov_z[-c(1:est_K0_plus_K2)] )
  est_K2 <- est_K0_plus_K2 - est_K0
  if(is.null(num_K2) == TRUE){
    if(est_K2 <= 0){
      cat("Estimation results in K0 >= K0+K2. K2 is set to 3.\n ")
      est_K2 <- 3
    }
  } else {
    est_K2 <- num_K2
  }

  # Estimation for v_k1 and lambda_k2
  est_lambda2 <- vector(mode="numeric", length = est_K2)

  V0 <- CoFA_CrossCov_fit$V
  P_V0 <- V0  %*% solve( t(V0) %*% V0  )%*% t(V0)
  P_V1 <- diag(1, p) - P_V0
  svd_indep_mul <- svd(P_V1 %*% (Cov_Z ) %*% P_V1)
  est_lambda2[1:est_K2] <- svd_indep_mul$d[1:est_K2] - est_error_omega_MLE



  ### 4. mBLUP
  ## 4-a. mMLUP for common components
  new_data <- rbind( t(as.matrix(y_sample) %*% W_sqrt ), t(z_sample_std) )
  mean_data <- rbind( t(y_centered %*% W_sqrt ), matrix(0, p, n) )
  x_demeaned <- ( new_data - mean_data  )

  cov_x_demeaned <- x_demeaned %*% t(x_demeaned) / (n-1)
  if(est_K0 ==1){
    M_trans <- sqrt( 1/(1+(est_beta[1])^2) ) * cbind( t( W_sqrt %*% B_tilde %*% CoFA_CrossCov_fit$U  ),  (est_beta[1]) * t(CoFA_CrossCov_fit$V)  )
    Blup_pred_com <- est_lambda0[1] * sqrt(1+ (est_beta[1])^2) * solve( M_trans %*% cov_x_demeaned %*% t(M_trans) ) %*% M_trans %*% x_demeaned
  } else {
    M_trans <- diag( sqrt( 1/( 1+(est_beta[1:est_K0])^2) ) ) %*% cbind( t( W_sqrt %*% B_tilde %*% CoFA_CrossCov_fit$U[, 1:est_K0]  ),  diag( (est_beta[1:est_K0]) ) %*% t(CoFA_CrossCov_fit$V[, 1:est_K0])  )
    Blup_pred_com <- diag( est_lambda0[1:est_K0]) %*% diag( sqrt( 1+ (est_beta[1:est_K0])^2 )) %*% solve( M_trans %*% cov_x_demeaned %*% t(M_trans) ) %*% M_trans %*% x_demeaned
  }



  ## 4-b. mMLUP for independent components in functional data
  y_demeaned_blup <- t(W_sqrt) %*% ( t(as.matrix(y_sample)) - t(y_centered)  )
  Phi_1 <- W_sqrt %*% B_tilde %*% svd_indep_fun$u[,1:est_K1]
  Blup_pred_ind_func <- diag( est_lambda1[1:est_K1] )   %*% solve( t(Phi_1) %*% y_demeaned_blup %*% t(y_demeaned_blup) %*% Phi_1 / (n-1) ) %*% t(Phi_1) %*% y_demeaned_blup


  ## 4-c. mMLUP for independent components in multivariate data
  Phi_2 <- svd_indep_mul$u[,1:est_K2]
  Blup_pred_ind_multi <- diag( est_lambda2[1:est_K2] )   %*% solve(t(Phi_2) %*% Cov_Z %*% Phi_2 ) %*% t(Phi_2) %*% t(z_sample_std)



  ### 5. Output
  phi_k0 <- B_tilde %*% CoFA_CrossCov_fit$U
  phi_k1 <- B_tilde %*% svd_indep_fun$u[,1:est_K1]
  v_k0 <- CoFA_CrossCov_fit$V
  v_k1 <- svd_indep_mul$u[,1:est_K2]

  result <- list(K0 = est_K0, K1 = est_K1, K2 = est_K2,
                 phi_k0 = phi_k0, phi_k1 = phi_k1, v_k0 = v_k0, v_k1 = v_k1,
                 lambda0 = est_lambda0, lambda1 = est_lambda1, lambda2 = est_lambda2,
                 beta = est_beta, sigma2_fun = est_error_epsilon_MLE, sigma2_mul = est_error_omega_MLE,
                 mBLUP_com = t(Blup_pred_com), mBLUP_fun = t(Blup_pred_ind_func), mBLUP_mul = t(Blup_pred_ind_multi),
                 tau_final = CoFA_CrossCov_fit$tau_final
                 )

  return(result)

}





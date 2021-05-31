#' A function to conduct real data application
#' Note that there is a seed number implanted for reproducibility
#' @export data_application
#' @param which_ear Either "right" or "left" to choose which ear of data is loaded
#' @param num_knot The vector of total numbers of knots
#' @param num_resampling_CV The number of resampling for the cross validation process with respect to tau
#' @return A list that contains results of estimation and prediction for different number of knots and the standardized multivariate data matrix
#' @importFrom face select.knots pspline
#' @importFrom splines splineDesign
#' @importFrom stats cov
#' @importFrom RMTstat qmp ptw
<<<<<<< HEAD
#' @import CVXR
=======
>>>>>>> 8641356d47df019f78f7131f416d14d7edfa9bd8


data_application <- function(which_ear = "right", num_knot = 12:44, num_resampling_CV = 100
){

  ### 1. Read in the real data
  real_data <- read_real_data(select_ear = which_ear)

  data_ready <- real_data$data_ready
  n <- dim(data_ready)[1]


  y_sample <- data_ready[, which(names(data_ready) == "WBXA001") : which(names(data_ready) == "WBXA107") ]
  t_frequencies <- data_ready[1, which(names(data_ready) == "WBXF001") : which(names(data_ready) == "WBXF107") ]

  t_temp  <- as.numeric(t_frequencies)
  t_temp <- log(t_temp)
  t_temp_std <- ( t_temp / (t_temp[length(t_temp)] - t_temp[1]) ) -  ( t_temp / (t_temp[length(t_temp)] - t_temp[1]) )[1]
  t_vec <- t_temp_std

  m <- length(t_vec)


  z_sample <- data_ready[, -c(1:323, dim(data_ready)[2]-1,  dim(data_ready)[2]) ]
  z_sample <- as.matrix(z_sample)
  z_sample_std <- scale(z_sample)
  z_sample_std_mat <- as.matrix(z_sample_std)
  p <- dim(z_sample_std)[2]

  # Two tuning parameters recommended in BEMA method
  BEMA_alpha <- 0.2
  BEMA_beta <- 0.1


  ### 2. Initialize objects to save estimation results
  iter_n <- length( num_knot )

  results_list <- vector(mode = "list", length = iter_n)

  est_K0 <- matrix(0,iter_n)
  est_K1 <- matrix(0,iter_n)
  est_K2 <- matrix(0,iter_n)

  est_eta <- matrix(0, iter_n, p)
  est_beta <- matrix(0, iter_n, p)
  est_lambda0 <- matrix(0, iter_n, p)
  est_lambda1 <- matrix(0, iter_n, m)
  est_lambda2 <- matrix(0, iter_n, p)

  est_error_epsilon_MLE <- matrix(0,iter_n,1)
  est_error_epsilon_BEMA <- matrix(0, iter_n, 1)
  est_error_epsilon_FACE <- matrix(0, iter_n, 1)
  est_error_epsilon_fpca_sc <- matrix(0, iter_n, 1)

  est_error_omega_MLE <- matrix(0, iter_n, 1)
  est_error_omega_BEMA <- matrix(0, iter_n, 1)

  set.seed(12345)

  ### 3. Start for loop for each number of knots
  for(i in 1:iter_n){

    cat(i, "th iteration begun \n")

    ### 4. Compute B_tilde matrix
    knots_set <- (num_knot - 7)[i]
    knots_vec <- select.knots(t_vec, knots = knots_set ,p=3)
    B_mat <- splineDesign(knots_vec, t_vec, ord = 4)

    # G <- (t(B_mat) %*% B_mat / m)
    w <- quadWeights(argvals = t_vec , method = "trapezoidal")
    W_mat <- diag(w)
    W_sqrt <- diag( sqrt(w) )
    G <- t(B_mat) %*% W_mat %*% B_mat
    G_sqrt <- eigen(G)$vectors %*% diag( sqrt(eigen(G)$values)  ) %*% t( eigen(G)$vectors )
    G_sqrt_inverse <- eigen(G)$vectors %*% diag( 1/sqrt(eigen(G)$values)  ) %*% t( eigen(G)$vectors )
    B_tilde <- B_mat %*% G_sqrt_inverse



    ### 5. Mean function estimation using P-spline
    T <- length(t_vec)
    id <- rep(1:n,each=T)
    t <- rep(t_vec,times=n)
    y_stack <- as.vector(t(y_sample))
    data <- data.frame(y=y_stack,
                       argvals = t,
                       subj = id)
    fit_pspline <- pspline(data,knots=knots_set)
    y_centered <- t(matrix(fit_pspline$fitted.values, length(t_vec), n))
    y_demeaned <- y_sample - y_centered
    y_demeaned_mat <- as.matrix(y_demeaned)





    ### 6. Estimation for sigma^{2}_{epsilon} and functional components from FACE method for comparison
    fit_face <- fpca_face_modified(Y=y_demeaned_mat,ydata=NULL,Y.pred = NULL,argvals=t_vec,pve = 0.99, npc  = NULL,
                                   var = TRUE, simul = FALSE, sim.alpha = 0.95,
                                   center=FALSE,knots_para=knots_set,p=3,m=2,lambda=NULL,alpha = 1,
                                   search.grid=TRUE,search.length=100,
                                   method="L-BFGS-B", lower=-20,upper=20, control=NULL)
    est_error_epsilon_FACE[i,1] <- fit_face$sigma2
    cov_y_face <- fit_face$Cov

    ### 4-2. fpca.sc
    fpca_sc_fit <- fpca_sc_modified(Y = y_demeaned_mat, ydata = NULL, Y.pred = NULL, argvals = t_vec, random.int = FALSE,
                                    nbasis = dim(B_tilde)[2] , pve = 0.99, npc = NULL, var = TRUE, simul = FALSE, sim.alpha = 0.95,
                                    useSymm = FALSE, makePD = FALSE, center = FALSE, cov.est.method = 2, integration = "trapezoidal")

    cov_fpcaSC_cov_hat <- fpca_sc_fit$cov.hat

    est_error_epsilon_fpca_sc[i,1] <- fpca_sc_fit$sigma2





    ### 7. CoFA for common component part
    CoFA_CrossCov_fit <- CoFA_CrossCov(Y_cent = y_demeaned_mat, Z_cent = z_sample_std_mat, B_tilde = B_tilde, t_vec=t_vec, gamma = 2,
                              control = list(num_rep= num_resampling_CV, num_fold=5, num_tau=100))

    Cross_cov <- 1/(n-1) *( t(y_demeaned_mat) %*% (z_sample_std_mat) )
    Gamma_mat <- t(B_tilde) %*% W_mat %*% Cross_cov

    svd_Gamma <- svd(Gamma_mat)


    K0_est <- CoFA_CrossCov_fit$rank
    est_K0[i] <- K0_est



    for(j in 1:K0_est){
      est_eta[i, j] <- t(CoFA_CrossCov_fit$U[,j]) %*% Gamma_mat %*% CoFA_CrossCov_fit$V[,j]
      est_beta[i,j] <- est_eta[i, j] / ( t(B_tilde %*% CoFA_CrossCov_fit$U[,j]) %*%  W_mat %*% cov_fpcaSC_cov_hat %*%  W_mat %*% B_tilde %*% CoFA_CrossCov_fit$U[,j] )
      est_lambda0[i, j] <-   (t(B_tilde %*% CoFA_CrossCov_fit$U[,j]) %*%  W_mat %*% cov_fpcaSC_cov_hat %*%  W_mat %*% B_tilde %*% CoFA_CrossCov_fit$U[,j]  )

    }





    ### 8. Estimation for independent component part
    ## 8-a functional data

    U0 <- CoFA_CrossCov_fit$U[,1:K0_est]
    P_U1 <- diag(1, dim(B_mat)[2]) - U0 %*% solve( t(U0) %*% U0) %*% t(U0)
    cov_y <- t(y_demeaned_mat) %*% (y_demeaned_mat) / (n-1)
    svd_indep_fun <- svd( P_U1 %*% t(B_tilde) %*% ( W_mat %*% cov_y %*% W_mat  ) %*% B_tilde  %*% t(P_U1) )


    for(j in 1:length(svd_indep_fun$d)){
      est_lambda1[i, j] <-   (t(B_tilde %*% svd_indep_fun$u[,j]) %*% W_mat %*% cov_fpcaSC_cov_hat %*% W_mat %*% B_tilde %*% svd_indep_fun$u[,j]  )
    }


    # Estimatiion K_{1} by BEMA and MLE of sigma^{2}_{epsilon}
    svd_Cov_y <- svd(cov_y)
    d_Cov_y <- svd_Cov_y$d

    p_tilde_fnt <- min(n, m)
    gamma_n_fnt <- m/n
    bulk_eigenvalue_order_fnt <- seq( ceiling(p_tilde_fnt * BEMA_alpha), floor( p_tilde_fnt * (1-BEMA_alpha) ) )
    bulk_eigenvalue_fnt <- d_Cov_y[bulk_eigenvalue_order_fnt]
    q_k_vec_fnt <-  qmp( bulk_eigenvalue_order_fnt/p_tilde_fnt , ndf=n, pdim=m, var=1, lower.tail = FALSE, log.p = FALSE )

    est_error_epsilon_BEMA[i,1] <- sum( bulk_eigenvalue_fnt * q_k_vec_fnt ) / sum(q_k_vec_fnt^(2))

    t_1minusbeta <- ptw( 1 - BEMA_beta , beta=1, lower.tail = TRUE, log.p = FALSE)
    T_hat_fnt <- est_error_epsilon_BEMA[i,1] *
      ( (1+sqrt(gamma_n_fnt))^(2) + t_1minusbeta * n ^(-2/3) * (gamma_n_fnt)^(-1/6) * (1+sqrt(gamma_n_fnt))^(4/3) )

    est_K0_plus_K1 <- sum( d_Cov_y > T_hat_fnt )
    est_error_epsilon_MLE[i,1]  <-   mean( d_Cov_y[-c(1:est_K0_plus_K1)] )
    est_K1[i] <- est_K0_plus_K1 - K0_est


    ## 8-b multivariate data
    # Estimatiion K_{2} by BEMA and MLE of sigma^{2}_{omega}
    Cov_Z <- 1/(n-1) * t(z_sample_std_mat) %*% z_sample_std_mat
    svd_Cov_z <- svd(Cov_Z)
    d_Cov_z <- svd_Cov_z$d

    p_tilde <- min(n, p)
    gamma_n <- p/n
    bulk_eigenvalue_order <- seq( ceiling(p_tilde * BEMA_alpha), floor( p_tilde * (1-BEMA_alpha) ) )
    bulk_eigenvalue <- d_Cov_z[bulk_eigenvalue_order]
    q_k_vec <-  qmp( bulk_eigenvalue_order/p_tilde , ndf=n, pdim=p, var=1, lower.tail = FALSE, log.p = FALSE )

    est_error_omega_BEMA[i,1] <- sum( bulk_eigenvalue * q_k_vec ) / sum(q_k_vec^(2))

    t_1minusbeta <- ptw( 1 - BEMA_beta , beta=1, lower.tail = TRUE, log.p = FALSE)
    T_hat <- est_error_omega_BEMA[i,1] *
      ( (1+sqrt(gamma_n))^(2) + t_1minusbeta * n ^(-2/3) * (gamma_n)^(-1/6) * (1+sqrt(gamma_n))^(4/3) )

    est_K0_plus_K2 <- sum( d_Cov_z > T_hat )
    est_error_omega_MLE[i,1]  <-   mean( d_Cov_z[-c(1:est_K0_plus_K2)] )
    est_K2[i] <- est_K0_plus_K2 - K0_est

    # Estimation for V1 and lambda_k2
    V0 <- CoFA_CrossCov_fit$V
    P_V0 <- V0  %*% solve( t(V0) %*% V0  )%*% t(V0)
    P_V1 <- diag(1, p) - P_V0
    svd_indep_mul <- svd(P_V1 %*% (Cov_Z ) %*% P_V1)
    est_lambda2[i, ] <- svd_indep_mul$d - est_error_omega_MLE[i,1]



    ### 9. mBLUP
    ## 9-a mMLUP for common components
    new_data <- rbind( t(as.matrix(y_sample) %*% W_sqrt ), t(z_sample_std_mat) )
    mean_data <- rbind( t(y_centered %*% W_sqrt ), matrix(0, p, n) )
    x_demeaned <- ( new_data - mean_data  )

    cov_x_demeaned <- x_demeaned %*% t(x_demeaned) / (n-1)
    if(K0_est ==1){
      M_trans <- sqrt( 1/(1+(est_beta[i,1])^2) ) * cbind( t( W_sqrt %*% B_tilde %*% CoFA_CrossCov_fit$U  ),  (est_beta[i,1]) * t(CoFA_CrossCov_fit$V)  )
      Blup_pred_com <- est_lambda0[i,1] * sqrt(1+ (est_beta[i, 1:K0_est])^2) * solve( M_trans %*% cov_x_demeaned %*% t(M_trans) ) %*% M_trans %*% x_demeaned
    } else {
      M_trans <- diag( sqrt( 1/( 1+(est_beta[i,1:K0_est])^2) ) ) %*% cbind( t( W_sqrt %*% B_tilde %*% CoFA_CrossCov_fit$U[, 1:K0_est]  ),  diag( (est_beta[i,1:K0_est]) ) %*% t(CoFA_CrossCov_fit$V[, 1:K0_est])  )
      Blup_pred_com <- diag( est_lambda0[i,1:K0_est]) %*% diag( sqrt( 1+ (est_beta[i, 1:K0_est])^2 )) %*% solve( M_trans %*% cov_x_demeaned %*% t(M_trans) ) %*% M_trans %*% x_demeaned
    }



    ## 9-b mMLUP for independent components in functional data
    y_demeaned_blup <- t(W_sqrt) %*% ( t(as.matrix(y_sample)) - t(y_centered)  )
    Phi_1 <- W_sqrt %*% B_tilde %*% svd_indep_fun$u[,1:3]
    Blup_pred_ind_func <- diag( est_lambda1[i,1:3] )   %*% solve( t(Phi_1) %*% y_demeaned_blup %*% t(y_demeaned_blup) %*% Phi_1 / (n-1) ) %*% t(Phi_1) %*% y_demeaned_blup


    ## 9-c mMLUP for independent components in multivariate data
    Phi_2 <- svd_indep_mul$u[,1:3]
    Blup_pred_ind_multi <- diag( est_lambda2[i,1:3] )   %*% solve(t(Phi_2) %*% Cov_Z %*% Phi_2 ) %*% t(Phi_2) %*% t(z_sample_std_mat)


    ### 10. Save all results
    results_list[[i]] <- list( "CoFA_CrossCov_fit"=CoFA_CrossCov_fit,
                               "svd_indep_fun"=svd_indep_fun, "svd_indep_mul"=svd_indep_mul,
                               "cov_y_face"=cov_y_face, "Cov_Z"=Cov_Z,
                               "fit_face"=fit_face,
                               "y_centered"= y_centered,
                               "Blup_pred_com"=Blup_pred_com,
                               "Blup_pred_ind_func"=Blup_pred_ind_func,
                               "Blup_pred_ind_multi"=Blup_pred_ind_multi,
                               "B_mat"=B_mat, "B_tilde"= B_tilde)
  }

  if(which_ear == "right"){
    save(results_list, y_sample, z_sample_std_mat, t_vec, t_frequencies,
         est_eta, est_beta, est_lambda0, est_lambda1, est_lambda2,
         est_K0, est_K1, est_K2,
         est_error_epsilon_fpca_sc, est_error_epsilon_BEMA, est_error_epsilon_FACE,est_error_epsilon_MLE,
         est_error_omega_MLE,
         est_error_omega_BEMA,
         file="real_data_app_right.RData")

  } else if(which_ear == "left"){
    save(results_list, y_sample, z_sample_std_mat, t_vec, t_frequencies,
         est_eta, est_beta, est_lambda0, est_lambda1, est_lambda2,
         est_K0, est_K1, est_K2,
         est_error_epsilon_fpca_sc, est_error_epsilon_BEMA, est_error_epsilon_FACE,est_error_epsilon_MLE,
         est_error_omega_MLE,
         est_error_omega_BEMA,
         file="real_data_app_left.RData")

  }

  res <- list("results_list"=results_list, "y_sample"=y_sample, "z_sample_std_mat"=z_sample_std_mat,
              "t_vec"=t_vec, "t_frequencies"=t_frequencies,
              "est_eta"=est_eta, "est_beta"=est_beta,
              "est_lambda0"=est_lambda0, "est_lambda1"=est_lambda1, "est_lambda2"=est_lambda2,
              "est_K0"=est_K0, "est_K1"=est_K1, "est_K2"=est_K2,
              "est_error_epsilon"=est_error_epsilon_MLE,
              "est_error_epsilon_BEMA"=est_error_epsilon_BEMA,
              "est_error_epsilon_FACE"=est_error_epsilon_FACE,
              "est_error_epsilon_fpca_sc"=est_error_epsilon_fpca_sc,
              "est_error_omega_MLE"=est_error_omega_MLE,
              "est_error_omega_BEMA"=est_error_omega_BEMA
  )

  return(res)
}


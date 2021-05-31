#' A function to conduct simulation study at a given signal noise level and a number of iterations
#' Note that two seed numbers are implanted in this function: the one for the true eigen vectors, and the other for reproducibility of data generation
#' @export simulation_study
#' @param snr The signal noise ratio value
#' @param sample_size A numeric vector of sample sizes
#' @param N The number of iterations
#' @param num_resampling_CV The number of repetitions for CoFA function
#' @return A list that contains 3 lists as elements for 3 sample sizes (100, 200, 500): each of them has N+1 elements for the first N having CoFA fits and BLUPS at a single iteration and the last one including the results of parameter estimation and assessment
#' @importFrom stats integrate cov rnorm quantile knots
#' @importFrom pracma randortho
#' @importFrom face select.knots pspline face.sparse
#' @importFrom splines splineDesign spline.des
#' @importFrom MASS mvrnorm
#' @importFrom Matrix Matrix
#' @importFrom RMTstat qmp ptw
<<<<<<< HEAD
#' @importFrom refund fpca.sc
#' @import CVXR
#' @import fdapace


simulation_study <- function(snr = 5, sample_size = c(100, 200, 500), N = 1000, num_resampling_CV = 100
){
=======


>>>>>>> 8641356d47df019f78f7131f416d14d7edfa9bd8

simulation_study <- function(snr = 5, sample_size = c(100, 200, 500), N = 1000, num_resampling_CV = 100
){

<<<<<<< HEAD
  results_final <- vector(mode = "list", length = 3)
  # A list of all reulsts of every scenario

  ###############################
  ###### Parameter Setting ######
  ###############################
  t_vec <- seq(from=0, to=1, by= 1/49)
  # observation points of functional data

  m <- length(t_vec)
  knots_set <- 16
  knots_vec <- select.knots(t_vec, knots = knots_set ,p=3)
  B_mat <- splineDesign(knots_vec, t_vec, ord = 4)
  # knot and B-spline basis setting
=======

  results_final <- vector(mode = "list", length = 3)
  # A list of all reulsts of every scenario

  ###############################
  ###### Parameter Setting ######
  ###############################
  t_vec <- seq(from=0, to=1, by= 1/49)
  # observation points of functional data
>>>>>>> 8641356d47df019f78f7131f416d14d7edfa9bd8

  m <- length(t_vec)
  knots_set <- 16
  knots_vec <- select.knots(t_vec, knots = knots_set ,p=3)
  B_mat <- splineDesign(knots_vec, t_vec, ord = 4)
  # knot and B-spline basis setting

  K0 <- 2
  # the number of common latent factors
  K1 <- 2
  # the number of independent latent factors in both functional and multivariate data

<<<<<<< HEAD
  lambdas <- 1/ (1:4 + 1)^2
  # variances of random factors for common and indepedent components in functional data

  sigma2_error <- sum(lambdas) / snr
  # variance of noise in functional data

  p <- 20
  # Dimension of multivariate data

  beta_k <- sqrt(c(3,2))
  # scale parameters in multivariate data

  length_sample_size <- length(sample_size)

  set.seed(12345)
  V_sqr <- randortho(p, type="orthonormal")
  # random orthonormal vectors for true eigen vectors in multivariate data
  V_mat <- V_sqr[,(1:K0)]
  V_rest <- V_sqr[, (K0+1):(K0+2)]
  # V_mat for common components and V_rest for independent ones

  Sigma_omega <- diag( c(0.4, 0.1))
  # variances of independent components in multivariate data

  error_omega <- 1/(snr * p) * ( sum( beta_k^2 * lambdas[c(1,3)]  ) +  sum( diag(Sigma_omega)) )
  # variance of random noise in multivariate data

  # Two tuning parameters recommended in BEMA method
  BEMA_alpha <- 0.2
  BEMA_beta <- 0.1
=======
  K0 <- 2
  # the number of common latent factors
  K1 <- 2
  # the number of independent latent factors in both functional and multivariate data

  lambdas <- 1/ (1:4 + 1)^2
  # variances of random factors for common and indepedent components in functional data

  sigma2_error <- sum(lambdas) / snr
  # variance of noise in functional data

  p <- 20
  # Dimension of multivariate data

  beta_k <- sqrt(c(3,2))
  # scale parameters in multivariate data

  length_sample_size <- length(sample_size)

  set.seed(12345)
  V_sqr <- randortho(p, type="orthonormal")
  # random orthonormal vectors for true eigen vectors in multivariate data
  V_mat <- V_sqr[,(1:K0)]
  V_rest <- V_sqr[, (K0+1):(K0+2)]
  # V_mat for common components and V_rest for independent ones

  Sigma_omega <- diag( c(0.4, 0.1))
  # variances of independent components in multivariate data

  error_omega <- 1/(snr * p) * ( sum( beta_k^2 * lambdas[c(1,3)]  ) +  sum( diag(Sigma_omega)) )
  # variance of random noise in multivariate data
>>>>>>> 8641356d47df019f78f7131f416d14d7edfa9bd8

  # Two tuning parameters recommended in BEMA method
  BEMA_alpha <- 0.2
  BEMA_beta <- 0.1

<<<<<<< HEAD
  ####################################################
  #### B_tilde orthonomal B-spline basis matrix  #####
  ####################################################
  # G <- (t(B_mat) %*% B_mat / m)
  w <- quadWeights(argvals = t_vec , method = "trapezoidal")
  W_mat <- diag(w)
  W_sqrt <- diag( sqrt(w) )
  G <- t(B_mat) %*% W_mat %*% B_mat
  G_sqrt <- eigen(G)$vectors %*% diag( sqrt(eigen(G)$values)  ) %*% t( eigen(G)$vectors )
  G_sqrt_inverse <- eigen(G)$vectors %*% diag( 1/sqrt(eigen(G)$values)  ) %*% t( eigen(G)$vectors )
  B_tilde <- B_mat %*% G_sqrt_inverse


  ###############################################
  ## local functions for simulation study only ##
  ###############################################

  #### Functions for data generation #####
  y <- function(i, t_y){
    mu(t_y) +
      t(sapply(1:K0,phi_k0,t_phi0=t_y)) %*% xi_sample[i, seq(from=1, to=K0+K1, by=2)] +
      t(sapply(1:K1,phi_k1,t_phi1=t_y)) %*% xi_sample[i, seq(from=2, to=K0+K1, by=2)]
  }
  # functional data generation without error

  z <- function(i){
    V_mat %*% (  beta_k * xi_sample[i, seq(from=1, to=K0+K1, by=2)]  )  + V_rest %*% xi_sample[i, seq(from=K0+K1+1, to=dim(xi_sample)[2], by=1)]
  }
  # multivariate data generation


=======

  ####################################################
  #### B_tilde orthonomal B-spline basis matrix  #####
  ####################################################
  # G <- (t(B_mat) %*% B_mat / m)
  w <- quadWeights(argvals = t_vec , method = "trapezoidal")
  W_mat <- diag(w)
  W_sqrt <- diag( sqrt(w) )
  G <- t(B_mat) %*% W_mat %*% B_mat
  G_sqrt <- eigen(G)$vectors %*% diag( sqrt(eigen(G)$values)  ) %*% t( eigen(G)$vectors )
  G_sqrt_inverse <- eigen(G)$vectors %*% diag( 1/sqrt(eigen(G)$values)  ) %*% t( eigen(G)$vectors )
  B_tilde <- B_mat %*% G_sqrt_inverse


  ###############################################
  ## local functions for simulation study only ##
  ###############################################

  #### Functions for data generation #####
  y <- function(i, t_y){
    mu(t_y) +
      t(sapply(1:K0,phi_k0,t_phi0=t_y)) %*% xi_sample[i, seq(from=1, to=K0+K1, by=2)] +
      t(sapply(1:K1,phi_k1,t_phi1=t_y)) %*% xi_sample[i, seq(from=2, to=K0+K1, by=2)]
  }
  # functional data generation without error

  z <- function(i){
    V_mat %*% (  beta_k * xi_sample[i, seq(from=1, to=K0+K1, by=2)]  )  + V_rest %*% xi_sample[i, seq(from=K0+K1+1, to=dim(xi_sample)[2], by=1)]
  }
  # multivariate data generation


>>>>>>> 8641356d47df019f78f7131f416d14d7edfa9bd8
  #### Functions for estimation performance assessment #####
  ## Cross covariance function
  est_cross_cov_estK0 <- function(cross_s,c_prime=K0_est, z_p=p_loop){
    est <- 0

    if(K0_est == 1){
      est <- as.numeric( est_eta[l,1] * ( t(CoFA_CrossCov_fit$U) %*%  G_sqrt_inverse %*%  t( splineDesign(knots_vec, cross_s , ord = 4)  ) ) ) * as.matrix(CoFA_CrossCov_fit$V)[z_p,1]
    } else {

      for(k in 1:c_prime){
        add <- as.numeric( est_eta[l,k] * ( t(CoFA_CrossCov_fit$U[,k]) %*%  G_sqrt_inverse %*%  t( splineDesign(knots_vec, cross_s , ord = 4)  ) ) ) * as.matrix(CoFA_CrossCov_fit$V)[z_p,k]
        est <- est+add
      }
    }
    return(est)
  }

  true_cross_cov <- function(cross_s,c_true=K0, z_p=p_loop){
    est <- 0
    for(k in 1:c_true){
      add <- as.numeric( beta_k[k] * lambdas[1+2*(k-1)] * phi_k0(k, cross_s) ) * V_mat[z_p,k]
      est <- est+add
    }
    return(est)
  }


  ## Eigen function
  phi_est_minus <- function(s){

    if(results_list[[i]]$CoFA_CrossCov_fit$rank==1){
      sqrd_minus <-  ( phi_k0(k, s) -
                         as.numeric( t(results_list[[i]]$CoFA_CrossCov_fit$U) %*%  G_sqrt_inverse %*%  t( splineDesign(knots_vec, s , ord = 4)  ) ) )^2

    } else {

      sqrd_minus <-  ( phi_k0(k, s) -
                         as.numeric( t(results_list[[i]]$CoFA_CrossCov_fit$U[,k]) %*%  G_sqrt_inverse %*%  t( splineDesign(knots_vec, s , ord = 4)  ) ) )^2
    }
    return(sqrd_minus)
  }


  phi_est_plus <- function(s){
    if(results_list[[i]]$CoFA_CrossCov_fit$rank==1){
      sqrd_plus <- ( phi_k0(k, s) +
                       as.numeric( t(results_list[[i]]$CoFA_CrossCov_fit$U) %*%  G_sqrt_inverse %*%  t( splineDesign(knots_vec, s , ord = 4)  ) ) )^2
    } else {

      sqrd_plus <- ( phi_k0(k, s) +
                       as.numeric( t(results_list[[i]]$CoFA_CrossCov_fit$U[,k]) %*%  G_sqrt_inverse %*%  t( splineDesign(knots_vec, s , ord = 4)  ) ) )^2
    }
    return(sqrd_plus)
  }

  phi1_est_minus <- function(s){
    sqrd_minus <-  ( phi_k1(k, s) -
                       as.numeric( t(results_list[[i]]$svd_indep_fun$u[,k]) %*%  G_sqrt_inverse %*%  t( splineDesign(knots_vec, s , ord = 4)  ) ) )^2
    return(sqrd_minus)
  }

  phi1_est_plus <- function(s){
    sqrd_plus <- ( phi_k1(k, s) +
                     as.numeric( t(results_list[[i]]$svd_indep_fun$u[,k]) %*%  G_sqrt_inverse %*%  t( splineDesign(knots_vec, s , ord = 4)  ) ) )^2
    return(sqrd_plus)
  }



  true_phi_10 <- rep(0,m)
  true_phi_20 <- rep(0,m)

  true_phi_11 <- rep(0,m)
  true_phi_21 <- rep(0,m)

  for(j in 1:m){
    true_phi_10[j] <- phi_k0(1, t_vec[j])
    true_phi_20[j] <- phi_k0(2, t_vec[j])

    true_phi_11[j] <- phi_k1(1, t_vec[j])
    true_phi_21[j] <- phi_k1(2, t_vec[j])
  }

  # true_cov_indep <- cbind(true_phi_11, true_phi_21) %*% diag(lambdas[c(2,4)]) %*% t(  cbind(true_phi_11, true_phi_21) )
  true_cov_y <- cbind(true_phi_10, true_phi_11, true_phi_20, true_phi_21) %*% diag(lambdas) %*% t(  cbind(true_phi_10, true_phi_11, true_phi_20, true_phi_21) )

  #############################
  ## Start simulation study ###
  #############################

  set.seed(123)
  for(s in 1:length_sample_size){

    ### Pre-assignment for record in the sth scenario ###
    results_list <- vector(mode = "list", length = N+1)
    # a list of results of the sth scenario
    n <- sample_size[s]
    # # of sample size


    est_eta <- matrix(0,N,p)
    est_beta <- matrix(0,N,p)
    est_lambda0 <- matrix(0,N,p)
    # variances of common components
    est_lambda1 <- matrix(0,N,m)
    # variances of independent components in functional data
    est_lambda2 <- matrix(0,N,m)
    # variances of independent components in multivariate data


    est_K0 <- matrix(0,N)
    est_K1 <- matrix(0,N)
    est_K2 <- matrix(0,N)

    est_error_epsilon_fpca_sc <- matrix(0, N, 1)
    est_error_epsilon_BEMA <- matrix(0, N, 1)
    est_error_epsilon_MLE <- matrix(0, N, 1)
    est_error_omega_BEMA <- matrix(0, N, 1)
    est_error_omega_MLE <- matrix(0, N, 1)

    cross_ise_estK0 <- matrix(0, N, p)
    # assessment of estimation performance for cross covariance


    ### Start For loop for a given sample size ###

    for(l in 1:N){

      cat(l,"th iteration begun \n")


      ### 1. Generating functional data
      total_num_components <- K0+K1+2
      zero_vec <- rep(0, total_num_components)

      xi_sample <- mvrnorm(n, mu=zero_vec, Sigma=diag(c(lambdas, diag(Sigma_omega) ) ) )
      error_mat_fnt <- matrix(rnorm(n * length(t_vec), 0, sd=sqrt(sigma2_error)), n, length(t_vec))
      error_mat_mult <- matrix(rnorm(n*p, 0,sd=sqrt(error_omega) ), n , p)

      y_sample_pre <- matrix(0,n,length(t_vec))

      for(i in 1:n){
        for(j in 1:length(t_vec)){
          y_sample_pre[i,j] <- y(i,t_vec[j])
        }
      }

      y_sample <- y_sample_pre + error_mat_fnt
      ## Graphical check for functional data sampling
      # plot(t_vec,mu(t_vec), type='l', lwd=3, col="red", xlim=c(0,1), ylim=c(-9,9), ylab="y(t)",xlab="t", main="Generated y(t), n=100")
      #  for(i in 1:n){
      #   lines(t_vec,results_final_1[[3]][[i]]$y_sample[i,])
      #  }
      # lines(t_vec,mu(t_vec), lwd=3, col="red")
      #  lines(t_vec, colMeans(y_sample), col="blue", lwd=3)


      ### 2. Generating multivariate data
      z_sample_pre <- matrix(0,n,p)
      for(i in 1:n){
        z_sample_pre[i,] <- z(i)
      }

      z_sample <- z_sample_pre + error_mat_mult


      ### 3. Estimation for mean function
      T <- length(t_vec)
      id <- rep(1:n,each=T)
      t <- rep(t_vec,times=n)
      y_stack <- as.vector(t(y_sample))
      data <- data.frame(y=y_stack,
                         argvals = t,
                         subj = id)
      fit_pspline <- pspline(data,knots=knots_set)

      y_centered <- t(matrix(fit_pspline$fitted.values, length(t_vec), n))
      # mean function matrix
      y_demeaned <- y_sample - y_centered
      # demeaned functional data




      ### 4. Covariance function and sigma^{2}_{epsilon} estimation by fpca.sc of refund package
      fpca_sc_fit <- fpca_sc_modified(Y = y_demeaned, ydata = NULL, Y.pred = NULL, argvals = t_vec, random.int = FALSE,
                                      nbasis = dim(B_tilde)[2] , pve = 0.99, npc = NULL, var = TRUE, simul = FALSE, sim.alpha = 0.95,
                                      useSymm = FALSE, makePD = FALSE, center = FALSE, cov.est.method = 2, integration = "trapezoidal")

      cov_fpcaSC_cov_hat <- fpca_sc_fit$cov.hat
      est_error_epsilon_fpca_sc[l,1] <- fpca_sc_fit$sigma2


      ### 5. Estimation for common component part
      CoFA_CrossCov_fit <- CoFA_CrossCov(Y_cent = y_demeaned, Z_cent = z_sample, B_tilde = B_tilde, t_vec=t_vec, gamma = 2,
                                control = list(num_rep=num_resampling_CV, num_fold=5, num_tau=100))

      K0_est <- CoFA_CrossCov_fit$rank
      est_K0[l] <- K0_est

      Cross_cov <- 1/(n-1) *( t(y_demeaned) %*% (z_sample) )
      Gamma_mat <- t(B_tilde) %*% W_mat %*% Cross_cov

      svd_Gamma <- svd(Gamma_mat)
      cov_y <- t(y_demeaned) %*% (y_demeaned) / (n-1)


      for(j in 1:K0_est){
        est_eta[l, j] <- t(CoFA_CrossCov_fit$U[,j]) %*% Gamma_mat %*% CoFA_CrossCov_fit$V[,j]
        est_beta[l, j] <- svd_Gamma$d[j] / (t(B_tilde %*% CoFA_CrossCov_fit$U[,j]) %*%  W_mat %*% cov_fpcaSC_cov_hat %*%  W_mat %*% B_tilde %*% CoFA_CrossCov_fit$U[,j]  )
        est_lambda0[l, j] <-   (t(B_tilde %*% CoFA_CrossCov_fit$U[,j]) %*%  W_mat %*% cov_fpcaSC_cov_hat %*%  W_mat %*% B_tilde %*% CoFA_CrossCov_fit$U[,j]  )
      }




      ### 6. Estimation for independent component part
      ## 6-a functional data

      U0 <- CoFA_CrossCov_fit$U[,1:K0_est]
      P_U1 <- diag(1, dim(B_mat)[2]) - U0 %*% solve( t(U0)%*% U0  )%*% t(U0)
      svd_indep_fun <- svd( P_U1 %*%  t(B_tilde) %*% ( W_mat %*% cov_y %*% W_mat  ) %*% B_tilde  %*% t(P_U1) )

      # The first 5 lambda_k1
      for(j in 1:5){
        est_lambda1[l, j] <-   (t(B_tilde %*% svd_indep_fun$u[,j]) %*% W_mat %*% cov_fpcaSC_cov_hat %*% W_mat %*% B_tilde %*% svd_indep_fun$u[,j]  )
      }

      # Estimatiion K_{1} by BEMA and MLE of sigma^{2}_{epsilon}
      cov_y <- t(y_demeaned) %*% (y_demeaned) / (n-1)
      svd_Cov_y <- svd(cov_y)
      d_Cov_y <- svd_Cov_y$d

      p_tilde_fnt <- min(n, m)
      gamma_n_fnt <- m/n
      bulk_eigenvalue_order_fnt <- seq( ceiling(p_tilde_fnt * BEMA_alpha), floor( p_tilde_fnt * (1-BEMA_alpha) ) )
      bulk_eigenvalue_fnt <- d_Cov_y[bulk_eigenvalue_order_fnt]
      q_k_vec_fnt <-  qmp( bulk_eigenvalue_order_fnt/p_tilde_fnt , ndf=n, pdim=m, var=1, lower.tail = FALSE, log.p = FALSE )

      est_error_epsilon_BEMA[l,1] <- sum( bulk_eigenvalue_fnt * q_k_vec_fnt ) / sum(q_k_vec_fnt^(2))
      t_1minusbeta <- ptw( 1 - BEMA_beta , beta=1, lower.tail = TRUE, log.p = FALSE)
      T_hat_fnt <- est_error_epsilon_BEMA[l,1] *
        ( (1+sqrt(gamma_n_fnt))^(2) + t_1minusbeta * n ^(-2/3) * (gamma_n_fnt)^(-1/6) * (1+sqrt(gamma_n_fnt))^(4/3) )

      est_K0_plus_K1 <- sum( d_Cov_y > T_hat_fnt )
      est_error_epsilon_MLE[l,1]  <-   mean( d_Cov_y[-c(1:est_K0_plus_K1)] )
      est_K1[l] <- est_K0_plus_K1 - K0_est


      ## 6-b multivariate data
      # Estimatiion K_{2} by BEMA and MLE of sigma^{2}_{omega}
      Cov_Z <- t(z_sample) %*% z_sample / (n-1)
      svd_Cov_z <- svd(Cov_Z)
      d_Cov_z <- svd_Cov_z$d

      p_tilde <- min(n, p)
      gamma_n <- p/n
      bulk_eigenvalue_order <- seq( ceiling(p_tilde * BEMA_alpha), floor( p_tilde * (1-BEMA_alpha) ) )
      bulk_eigenvalue <- d_Cov_z[bulk_eigenvalue_order]
      q_k_vec <-  qmp( bulk_eigenvalue_order/p_tilde , ndf=n, pdim=p, var=1, lower.tail = FALSE, log.p = FALSE )

      est_error_omega_BEMA[l,1] <- sum( bulk_eigenvalue * q_k_vec ) / sum(q_k_vec^(2))

      t_1minusbeta <- ptw( 1 - BEMA_beta , beta=1, lower.tail = TRUE, log.p = FALSE)
      T_hat <- est_error_omega_BEMA[l,1] *
        ( (1+sqrt(gamma_n))^(2) + t_1minusbeta * n ^(-2/3) * (gamma_n)^(-1/6) * (1+sqrt(gamma_n))^(4/3) )

      est_K0_plus_K2 <- sum( d_Cov_z > T_hat )
      est_error_omega_MLE[l,1]  <-   mean( d_Cov_z[-c(1:est_K0_plus_K2)] )
      est_K2[l] <- est_K0_plus_K2 - K0_est


      # Estimation for V1 and lambda_k2
      V0 <- CoFA_CrossCov_fit$V
      P_V0 <- V0  %*% solve( t(V0) %*% V0  )%*% t(V0)
      P_V1 <- diag(1, p) - P_V0

      svd_indep_mul <- svd(P_V1 %*% (Cov_Z ) %*% P_V1)

      est_lambda2[l, 1:(length(svd_indep_mul$d))] <- diag( t(svd_indep_mul$u) %*% Cov_Z %*% svd_indep_mul$u ) - est_error_omega_MLE[l,1]






      ### 7. Marginal and Joint BLUPs
      ## 7-a xi_k0 jBLUP
      if(K0_est ==1){

        Lambda0_product_N_transpose <- cbind(  est_lambda0[l, 1:K0_est] *  t( W_sqrt %*% B_tilde %*% CoFA_CrossCov_fit$U  ),  est_eta[l, 1:K0_est] *  t(CoFA_CrossCov_fit$V) )


      } else {

        Lambda0_product_N_transpose <- cbind(   diag(est_lambda0[l, 1:K0_est]) %*% t(W_sqrt %*% B_tilde %*% CoFA_CrossCov_fit$U[, 1:K0_est]  ),  diag(est_eta[l, 1:K0_est]) %*%  t(CoFA_CrossCov_fit$V[, 1:K0_est]) )

      }
      new_data <- rbind(  t(y_sample %*% W_sqrt ), t(z_sample) )
      mean_data <- rbind(  t(y_centered %*% W_sqrt ), matrix(0, p, n) )
      V_hat_x_sample <- (new_data - mean_data ) %*% t(new_data - mean_data ) / (n-1)
      jBlup_pred_x <- Lambda0_product_N_transpose %*% solve(V_hat_x_sample) %*% ( new_data - mean_data  )


      ## 7-b xi_k0  mBLUP
      new_data <- rbind( t(y_sample %*% W_sqrt ), t(z_sample) )
      mean_data <- rbind( t(y_centered %*% W_sqrt ), matrix(0, p, n) )
      x_demeaned <- ( new_data - mean_data  )
      cov_x_demeaned <- x_demeaned %*% t(x_demeaned) / (n-1)

      if(K0_est ==1){
        M_trans <- sqrt( 1/(1+(est_beta[l,1])^2) ) * cbind( t( W_sqrt %*% B_tilde %*% CoFA_CrossCov_fit$U  ),  (est_beta[l,1]) * t(CoFA_CrossCov_fit$V)  )
        mBlup_pred_x <- est_lambda0[l,1] * sqrt(1+ (est_beta[l, 1:K0_est])^2) * solve( M_trans %*% cov_x_demeaned %*% t(M_trans) ) %*% M_trans %*% x_demeaned
      } else {
        M_trans <- diag( sqrt( 1/( 1+(est_beta[l,1:K0_est])^2) ) ) %*% cbind( t( W_sqrt %*% B_tilde %*% CoFA_CrossCov_fit$U[, 1:K0_est]  ),  diag( (est_beta[l,1:K0_est]) ) %*% t(CoFA_CrossCov_fit$V[, 1:K0_est])  )
        mBlup_pred_x <- diag( est_lambda0[l,1:K0_est]) %*% diag( sqrt( 1+ (est_beta[l, 1:K0_est])^2 )) %*% solve( M_trans %*% cov_x_demeaned %*% t(M_trans) ) %*% M_trans %*% x_demeaned
      }


      ## 7-c xi_k1 jBLUP
      cov_y_sample <- t(y_demeaned %*% W_sqrt ) %*% (y_demeaned %*% W_sqrt ) / (n-1)
      Lambda0_product_Phi1_transpose <-   diag(est_lambda1[l, 1:2]) %*%  t( W_sqrt %*% B_tilde %*% svd_indep_fun$u[,1:2]  )
      jBlup_pred_y <- Lambda0_product_Phi1_transpose %*% solve(cov_y_sample) %*% ( ( t(y_sample %*% W_sqrt ) - t(y_centered %*% W_sqrt ) )  )


      ## 7-d xi_k1 mBLUP
      Phi_1 <- W_sqrt %*% B_tilde %*% svd_indep_fun$u[,1:2]
      mBlup_pred_y <- diag(est_lambda1[l, 1:2]) %*% solve( t(Phi_1) %*% cov_y_sample %*% Phi_1  ) %*% t(Phi_1) %*% ( t(W_sqrt) %*% ( t(y_sample) - t(y_centered) ) )


      ## 7-e xi_k2 jBLUP and mBLUP
      cov_z_sample <- t(z_sample) %*% (z_sample) / (n-1)
      V1_est <- svd_indep_mul$u[,1:2]

      mBlup_pred_z <- diag(est_lambda2[l,1:2]) %*% solve( t(V1_est) %*% cov_z_sample %*% V1_est ) %*% t(V1_est) %*% t(z_sample)
      jBlup_pred_z <- diag(est_lambda2[l,1:2]) %*% t(V1_est) %*% solve(cov_z_sample)  %*% t(z_sample)



      ### 8. Assessment of estimation performance for cross-covariance
      for(p_loop in 1:p){
        cross_cov_diff<- function(s){
          sqrd <- (est_cross_cov_estK0(cross_s=s) - true_cross_cov(cross_s=s))^2

          return(sqrd)
        }


        if(K0_est == 1){
          cross_ise <- integrate(cross_cov_diff, lower=0, upper=1)
          cross_ise_estK0[l,p_loop]<- cross_ise$value
        } else {
          cross_ise <- integrate(cross_cov_diff, lower=0, upper=1)
          cross_ise_estK0[l,p_loop]<- cross_ise$value
        }

      }






      ### 9. Save results at the lth iteration
      results_list[[l]] <- list(
        "xi_sample"=xi_sample,
        "Cov_Z"=Cov_Z,
        "cov_y"=cov_y,
        "svd_Gamma"=svd_Gamma,

        "CoFA_CrossCov_fit" = CoFA_CrossCov_fit,

        "svd_indep_fun"=svd_indep_fun,
        "svd_indep_mul"=svd_indep_mul,

        "mBlup_pred_x"= mBlup_pred_x,
        "jBlup_pred_x" = jBlup_pred_x,

        "mBlup_pred_y" = mBlup_pred_y,
        "jBlup_pred_y" = jBlup_pred_y,

        "mBlup_pred_z" = mBlup_pred_z,
        "jBlup_pred_z" = jBlup_pred_z

      )



    }
    # End of N iterations


    ### 10. Assessment of estimation performance for eigen functions and vectors
    ise_mat <- matrix(0,N,3)
    ise1_mat <- matrix(0,N,3)
    se_mat <- matrix(0,N,3)
    se1_mat <- matrix(0,N,3)


    for(i in 1:N){

      if(results_list[[i]]$CoFA_CrossCov_fit$rank == 1){
        k <- 1
        phi_minus <- integrate(phi_est_minus, lower=0, upper=1)
        phi_plus <- integrate(phi_est_plus, lower=0, upper=1)


        ise_mat[i,k]<- min(phi_minus$value, phi_plus$value)


        nu_est <- as.matrix(results_list[[i]]$CoFA_CrossCov_fit$V)

        diff_vec_minus <- sum( ( V_mat[,k] -   nu_est[,k]  )^2  )
        diff_vec_plus <- sum( ( V_mat[,k] +   nu_est[,k]  )^2  )

        se_mat[i, k] <- min(diff_vec_minus, diff_vec_plus)

      } else {
        for(k in 1:K0){

          phi_minus <- integrate(phi_est_minus, lower=0, upper=1)
          phi_plus <- integrate(phi_est_plus, lower=0, upper=1)


          ise_mat[i,k]<- min(phi_minus$value, phi_plus$value)

        }

        nu_est <- as.matrix(results_list[[i]]$CoFA_CrossCov_fit$V)
        for(k in 1:K0){


          diff_vec_minus <- sum( ( V_mat[,k] -   nu_est[,k]  )^2  )
          diff_vec_plus <- sum( ( V_mat[,k] +   nu_est[,k]  )^2  )

          se_mat[i, k] <- min(diff_vec_minus, diff_vec_plus)
        }

      } ## conditional if ends

      ## for phi_k1 and nu_k1
      for(k in 1:2){
        phi1_minus <- integrate(phi1_est_minus, lower=0, upper=1)
        phi1_plus <- integrate(phi1_est_plus, lower=0, upper=1)


        ise1_mat[i,k]<- min(phi1_minus$value, phi1_plus$value)

        nu_k1_est <- as.matrix(results_list[[i]]$svd_indep_mul$u)


        diff_vec_k1_minus <- sum( ( V_rest[,k] -   nu_k1_est[,k]  )^2  )
        diff_vec_k1_plus <- sum( ( V_rest[,k] +   nu_k1_est[,k]  )^2  )

        se1_mat[i, k] <- min(diff_vec_k1_minus, diff_vec_k1_plus)

      }
    } # 1 to N for loop end



    ### 11. Save all results of the sth scenario in the (N+1)th element

    results_list[[N+1]] <- list("est_K0"=est_K0,
                                "est_K1"=est_K1,
                                "est_K2"=est_K2,

                                "est_eta"=est_eta,
                                "est_beta"=est_beta,
                                "est_lambda0"=est_lambda0,
                                "est_lambda1"=est_lambda1,
                                "est_lambda2"=est_lambda2,

                                "ise_mat"=ise_mat,
                                "se_mat"=se_mat,
                                "ise1_mat"=ise1_mat,
                                "se1_mat"=se1_mat,
                                "cross_ise_estK0"=cross_ise_estK0,

                                "true_cov_y"=true_cov_y,
                                "true_phi_10"=true_phi_10,
                                "true_phi_20"=true_phi_20,
                                "true_phi_11"=true_phi_11,
                                "true_phi_21"=true_phi_21,


                                "est_error_epsilon_fpca_sc"=est_error_epsilon_fpca_sc,
                                "est_error_epsilon_BEMA"=est_error_epsilon_BEMA,
                                "est_error_epsilon_MLE"=est_error_epsilon_MLE,
                                "est_error_omega_BEMA"=est_error_omega_BEMA,
                                "est_error_omega_MLE"=est_error_omega_MLE

    )


    ### 12. Save all results of every scenario
    results_final[[s]] <- results_list
    rm(list=c("results_list"))
  }
  # end of all simulation
  save(results_final, file="sim_results.RData")
  return(results_final)
}

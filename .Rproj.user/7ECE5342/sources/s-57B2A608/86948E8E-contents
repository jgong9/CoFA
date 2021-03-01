#' A function to conduct simulation study at a given signal noise level and a number of iterations
#' Note that two seed numbers are implanted in this function: the one for the true eigen vectors, and the other for reproducibility of data generation
#' @export simulation_study
#' @param snr The signal noise ratio value
#' @param N The number of iterations
#' @return A list that contains 3 lists as elements for 3 sample sizes (100, 200, 500): each of them has N+1 elements for the first N having CoFA fits and BLUPS at a single iteration and the last one including the results of parameter estimation and assessment
#' @importFrom stats integrate cov rnorm quantile knots
#' @importFrom pracma randortho
#' @importFrom face select.knots pspline
#' @importFrom splines splineDesign spline.des
#' @importFrom MASS mvrnorm
#' @importFrom Matrix Matrix
#' @import CVXR


simulation_study <- function(snr = 5, N = 100){


results_final <- vector(mode = "list", length = 3)
# A list of all reulsts of every scenario

###### Parameter Setting #############
t_vec <- seq(from=0, to=1, by= 1/49)
# observation points of functional data

m <- length(t_vec)
knots_set <- 16
knots_vec <- select.knots(t_vec, knots = knots_set ,p=3)
B_mat <- splineDesign(knots_vec, t_vec, ord = 4)
# knot and B-spline basis setting



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

####################################################
#### B_tilde orthonomal B-spline basis matrix  #####
####################################################
G <- (t(B_mat) %*% B_mat / m)
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


#### Functions for estimation performance assessment #####
## Cross covariance function
est_cross_cov_estK0 <- function(cross_s,c_prime=K0_est, z_p=p_loop){
  est <- 0

  if(K0_est == 1){
    est <- as.numeric( CoFA_fit$D * ( t(CoFA_fit$U) %*%  G_sqrt_inverse %*%  t( splineDesign(knots_vec, cross_s , ord = 4)  ) ) ) * as.matrix(CoFA_fit$V)[z_p,1]
  } else {

    for(k in 1:c_prime){
      add <- as.numeric( diag(CoFA_fit$D)[k] * ( t(CoFA_fit$U[,k]) %*%  G_sqrt_inverse %*%  t( splineDesign(knots_vec, cross_s , ord = 4)  ) ) ) * as.matrix(CoFA_fit$V)[z_p,k]
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

  if(results_list[[i]]$CoFA_fit$rank==1){
    sqrd_minus <-  ( phi_k0(k, s) -
                       as.numeric( t(results_list[[i]]$CoFA_fit$U) %*%  G_sqrt_inverse %*%  t( splineDesign(knots_vec, s , ord = 4)  ) ) )^2

  } else {

    sqrd_minus <-  ( phi_k0(k, s) -
                       as.numeric( t(results_list[[i]]$CoFA_fit$U[,k]) %*%  G_sqrt_inverse %*%  t( splineDesign(knots_vec, s , ord = 4)  ) ) )^2
  }
  return(sqrd_minus)
}


phi_est_plus <- function(s){
  if(results_list[[i]]$CoFA_fit$rank==1){
    sqrd_plus <- ( phi_k0(k, s) +
                     as.numeric( t(results_list[[i]]$CoFA_fit$U) %*%  G_sqrt_inverse %*%  t( splineDesign(knots_vec, s , ord = 4)  ) ) )^2
  } else {

    sqrd_plus <- ( phi_k0(k, s) +
                     as.numeric( t(results_list[[i]]$CoFA_fit$U[,k]) %*%  G_sqrt_inverse %*%  t( splineDesign(knots_vec, s , ord = 4)  ) ) )^2
  }
  return(sqrd_plus)
}

phi1_est_minus <- function(s){
  sqrd_minus <-  ( phi_k1(k, s) -
                     as.numeric( t(results_list[[i]]$svd_k1$u[,k]) %*%  G_sqrt_inverse %*%  t( splineDesign(knots_vec, s , ord = 4)  ) ) )^2
  return(sqrd_minus)
}

phi1_est_plus <- function(s){
  sqrd_plus <- ( phi_k1(k, s) +
                   as.numeric( t(results_list[[i]]$svd_k1$u[,k]) %*%  G_sqrt_inverse %*%  t( splineDesign(knots_vec, s , ord = 4)  ) ) )^2
  return(sqrd_plus)
}





#############################
## Start simulation study ###
#############################

set.seed(123)
for(s in 1:3){

### Pre-assignment for record in the sth scenario ###
results_list <- vector(mode = "list", length = N+1)
# a list of results of the sth scenario
n <- c(100, 200, 500)[s]
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
est_error_epsilon <- matrix(0,N,1)
est_error_omega <- matrix(0, N, 1)

cross_ise_estK0 <- matrix(0, N, p)
# assessment of estimation performance for cross covariance
# cross_ise_estK0_2







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

  Cross_cov <- cov(y_demeaned, z_sample)
  Target_mat <- t(B_tilde) %*% Cross_cov / m

  ### 4. Estimation for sigma^{2}_{epsilon}
  fit_face <- fpca_face_adjust(Y=y_demeaned,ydata=NULL,Y.pred = NULL,argvals=t_vec, pve = 0.99, npc  = NULL,
                               var = TRUE, simul = FALSE, sim.alpha = 0.95,
                               center=FALSE,knots_para=knots_set,p=3,m=2,lambda=NULL,alpha = 1,
                               search.grid=TRUE,search.length=100,
                               method="L-BFGS-B", lower=-20,upper=20, control=NULL)
  est_error_epsilon[l,1] <- fit_face$sigma2



  ### 5. Estimation for common component part
  CoFA_fit <- CoFA_CrossCov(Y_cent = y_demeaned, Z_cent = z_sample, B_tilde = B_tilde, gamma = 2,
                            control = list(num_rep=100, num_fold=5, num_tau=100))

  cov_y <- t(y_demeaned) %*% (y_demeaned) / (n-1)
  C_tilde <- 1/(m^2) * t(B_tilde) %*% ( cov_y - diag(rep(est_error_epsilon[l,1],m)  ) ) %*% B_tilde

  K0_est <- CoFA_fit$rank
  est_K0[l] <- K0_est

  if(K0_est == 1){
    est_eta[l, 1] <- CoFA_fit$D
    est_beta[l,1] <- est_eta[l, 1] / (t(CoFA_fit$U) %*% C_tilde %*% CoFA_fit$U   )
    est_lambda0[l, 1] <-  (t(CoFA_fit$U) %*% C_tilde %*% CoFA_fit$U  )

  } else {
    est_eta[l,1:length(diag(CoFA_fit$D))] <- diag(CoFA_fit$D)
    for(j in 1:K0_est){
      est_beta[l,j] <- diag(CoFA_fit$D)[j] / (t(CoFA_fit$U[,j]) %*% C_tilde %*% CoFA_fit$U[,j]   )
      est_lambda0[l, j] <-  (t(CoFA_fit$U[,j]) %*%  C_tilde  %*% CoFA_fit$U[,j]  )
    }
  }



  ### 6. Estimation for independent component part
  ## 6-a functional data
  if(K0_est == 1){

    A0 <- CoFA_fit$U
    P1 <- diag(1, dim(B_mat)[2]) - A0 %*% solve( t(A0)%*% A0  )%*% t(A0)
    svd_k1 <- svd( P1 %*%  C_tilde  %*% t(P1) )

  } else{

    A0 <- CoFA_fit$U[,1:K0_est]
    P1 <- diag(1, dim(B_mat)[2]) - A0 %*% solve( t(A0)%*% A0  )%*% t(A0)
    svd_k1 <- svd( P1 %*%  C_tilde %*% t(P1) )

  }

  est_lambda1[l,1:length(svd_k1$d)] <- svd_k1$d


  ## 6-b multivariate data
  Cov_Z <- t(z_sample) %*% z_sample / (n-1)
  V0 <- CoFA_fit$V
  P_v0 <- V0  %*% solve( t(V0) %*% V0  )%*% t(V0)

  # estimation for sigma^{2}_{omega}
  # based on Ordinary Least Square with a constraint
  var_omega <- var_omega_est(z_data_cent = z_sample, V0_mat = CoFA_fit$V,
                                 est_beta_vec =  est_beta[l, 1:K0_est],
                                 est_eta_vec = est_eta[l, 1:K0_est],
                                 est_K0 = K0_est)

  est_error_omega[l,1]  <-   var_omega
  # end of sigma^{2}_{omega} estimation

  P_v1 <- diag(1, p) - P_v0
  svd_nu_k1 <- svd(P_v1 %*% (Cov_Z - diag(est_error_omega[l,1] , p)) %*% P_v1)
  est_lambda2[l, 1:(length(svd_nu_k1$d))] <- svd_nu_k1$d


  ###############################################################################


  ### 7. Marginal and Joint BLUPs
  ## 7-a xi_k0 jBLUP
  if(K0_est ==1){

    Lambda0_product_N_transpose <- cbind(  est_lambda0[l, 1:K0_est] * 1/sqrt(m) * t( B_tilde %*% CoFA_fit$U  ),  est_eta[l, 1:K0_est] *  t(CoFA_fit$V) )


  } else {

    Lambda0_product_N_transpose <- cbind(  1/sqrt(m) * diag(est_lambda0[l, 1:K0_est]) %*% t( B_tilde %*% CoFA_fit$U[, 1:K0_est]  ),  diag(est_eta[l, 1:K0_est]) %*%  t(CoFA_fit$V[, 1:K0_est]) )

  }
  new_data <- rbind( 1/sqrt(m) * t(y_sample), t(z_sample) )
  mean_data <- rbind( 1/sqrt(m) * t(y_centered), matrix(0, p, n) )
  V_hat_x_sample <- (new_data - mean_data ) %*% t(new_data - mean_data ) / (n-1)
  jBlup_pred_x <- Lambda0_product_N_transpose %*% solve(V_hat_x_sample) %*% ( new_data - mean_data  )


  ## 7-b xi_k0  mBLUP
  new_data <- rbind( 1/sqrt(m) * t(y_sample), t(z_sample) )
  mean_data <- rbind( 1/sqrt(m) * t(y_centered), matrix(0, p, n) )
  x_demeaned <- ( new_data - mean_data  )
  cov_x_demeaned <- x_demeaned %*% t(x_demeaned) / (n-1)

  if(K0_est ==1){
    M_mat <- sqrt( 1/(1+(est_beta[l,1])^2) ) * cbind( 1/sqrt(m) * t( B_tilde %*% CoFA_fit$U  ),  (est_beta[l,1]) * t(CoFA_fit$V)  )
    x_new <- M_mat %*% x_demeaned
    mBlup_pred_x <- est_lambda0[l,1] * sqrt(1+ (est_beta[l, 1:K0_est])^2) * solve( M_mat %*% cov_x_demeaned %*% t(M_mat) ) %*% x_new
  } else {
    M_mat <- diag( sqrt( 1/( 1+(est_beta[l,1:K0_est])^2) ) ) %*% cbind( 1/sqrt(m) * t( B_tilde %*% CoFA_fit$U[, 1:K0_est]  ),  diag( (est_beta[l,1:K0_est]) ) %*% t(CoFA_fit$V[, 1:K0_est])  )
    x_new <- M_mat %*% x_demeaned
    mBlup_pred_x <- diag( est_lambda0[l,1:K0_est]) %*% diag( sqrt( 1+ (est_beta[l, 1:K0_est])^2 )) %*% solve( M_mat %*% cov_x_demeaned %*% t(M_mat) ) %*% x_new
  }


  ## 7-c xi_k1 jBLUP
  cov_y_sample <- 1/m * t(y_demeaned) %*% (y_demeaned) / (n-1)
  Lambda0_product_Phi1_transpose <-   1/sqrt(m) * diag(est_lambda1[l, 1:2]) %*%  t( B_tilde %*% svd_k1$u[,1:2]  )
  jBlup_pred_y <- Lambda0_product_Phi1_transpose %*% solve(cov_y_sample) %*% ( 1/sqrt(m) * ( t(y_sample) - t(y_centered) )  )


  ## 7-d xi_k1 mBLUP
  Phi_1 <- 1/sqrt(m) * B_tilde %*% svd_k1$u[,1:2]
  mBlup_pred_y <- diag(est_lambda1[l, 1:2]) %*% solve( t(Phi_1) %*% cov_y_sample %*% Phi_1  ) %*% t(Phi_1) %*% ( 1/sqrt(m) * ( t(y_sample) - t(y_centered) ) )


  ## 7-e xi_k2 jBLUP and mBLUP
  cov_z_sample <- t(z_sample) %*% (z_sample) / (n-1)
  V_1_est <- svd_nu_k1$u[,1:2]
  est_lambda_2 <- svd_nu_k1$d[1:2]

  mBlup_pred_z <- diag(est_lambda_2) %*% solve( t(V_1_est) %*% cov_z_sample %*% V_1_est ) %*% t(V_1_est) %*% t(z_sample)
  jBlup_pred_z <- diag(est_lambda_2) %*% t(V_1_est) %*% solve(cov_z_sample)  %*% t(z_sample)



  ### 8. Assessment of estimation performance for cross-covariance
  for(p_loop in 1:p){
    cross_var_mise1<- function(s){
      sqrd <- (est_cross_cov_estK0(cross_s=s) - true_cross_cov(cross_s=s))^2

      return(sqrd)
    }


    if(K0_est == 1){
      cross_ise1 <- integrate(cross_var_mise1, lower=0, upper=1)
      cross_ise_estK0[l,p_loop]<- cross_ise1$value
    } else {
      cross_ise1 <- integrate(cross_var_mise1, lower=0, upper=1)
      cross_ise_estK0[l,p_loop]<- cross_ise1$value
    }

  }






  ### 9. Save results at the lth iteration
  results_list[[l]] <- list(
    "xi_sample"=xi_sample,

    "CoFA_fit" = CoFA_fit,
    "svd_k1"=svd_k1,
    "svd_nu_k1"=svd_nu_k1,
    "jBlup_pred_x" = jBlup_pred_x,

    "mBlup_pred_y" = mBlup_pred_y,
    "jBlup_pred_y" = jBlup_pred_y,

    "mBlup_pred_x"= mBlup_pred_x,

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

  if(results_list[[i]]$CoFA_fit$rank == 1){
    k <- 1
    phi_minus <- integrate(phi_est_minus, lower=0, upper=1)
    phi_plus <- integrate(phi_est_plus, lower=0, upper=1)


    ise_mat[i,k]<- min(phi_minus$value, phi_plus$value)


    nu_est <- as.matrix(results_list[[i]]$CoFA_fit$V)

    diff_vec_minus <- sum( ( V_mat[,k] -   nu_est[,k]  )^2  )
    diff_vec_plus <- sum( ( V_mat[,k] +   nu_est[,k]  )^2  )

    se_mat[i, k] <- min(diff_vec_minus, diff_vec_plus)

  } else {
    for(k in 1:K0){

      phi_minus <- integrate(phi_est_minus, lower=0, upper=1)
      phi_plus <- integrate(phi_est_plus, lower=0, upper=1)


      ise_mat[i,k]<- min(phi_minus$value, phi_plus$value)

    }

    nu_est <- as.matrix(results_list[[i]]$CoFA_fit$V)
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

    nu_k1_est <- as.matrix(results_list[[i]]$svd_nu_k1$u)


    diff_vec_k1_minus <- sum( ( V_rest[,k] -   nu_k1_est[,k]  )^2  )
    diff_vec_k1_plus <- sum( ( V_rest[,k] +   nu_k1_est[,k]  )^2  )

    se1_mat[i, k] <- min(diff_vec_k1_minus, diff_vec_k1_plus)


  }

} # 1 to N for loop end




### 11. Save all results of the sth scenario in the (N+1)th element

results_list[[N+1]] <- list("est_K0"=est_K0, "est_eta"=est_eta, "est_beta"=est_beta, "est_lambda0"=est_lambda0,
                            "est_lambda1"=est_lambda1, "ise_mat"=ise_mat, "se_mat"=se_mat, "ise1_mat"=ise1_mat, "se1_mat"=se1_mat,
                            "est_lambda2"=est_lambda2,
                            "est_error_epsilon"=est_error_epsilon, "est_error_omega"=est_error_omega,
                            "cross_ise_estK0"=cross_ise_estK0
)


### 12. Save all results of every scenario
results_final[[s]] <- results_list
rm(list=c("results_list"))
}
# end of all simulation
save(results_final, file="sim_results.RData")
return(results_final)
}

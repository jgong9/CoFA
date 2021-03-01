#' A function to conduct real data application
#' Note that there is a seed number implanted for reproducibility
#' @export data_application
#' @param which_ear Either "right" or "left" to choose which ear of data is loaded
#' @param num_knot The vector of total numbers of knots
#' @return A list that contains results of estimation and prediction for different number of knots and the standardized multivariate data matrix
#' @importFrom face select.knots pspline
#' @importFrom splines splineDesign
#' @importFrom stats cov
#' @import CVXR


data_application <- function(which_ear = "right", num_knot = 12:44){

### 1. Read in the real data
real_data <- read_real_data(select_ear = which_ear)

data_ready <- real_data$data_ready
n <- dim(data_ready)[1]


y_sample <- data_ready[, which(names(data_ready) == "WBXA001") : which(names(data_ready) == "WBXA107") ]
t_frequencies <- data_ready[1, which(names(data_ready) == "WBXF001") : which(names(data_ready) == "WBXF107") ]

t_temp  <- as.numeric(t_frequencies)
t_temp_std <- ( t_temp / (t_temp[length(t_temp)] - t_temp[1]) ) -  ( t_temp / (t_temp[length(t_temp)] - t_temp[1]) )[1]
t_vec <- t_temp_std
m <- length(t_vec)


z_sample <- data_ready[, -c(1:323, dim(data_ready)[2]-1,  dim(data_ready)[2]) ]
z_sample <- as.matrix(z_sample)
z_sample_std <- scale(z_sample)
z_sample_std_mat <- as.matrix(z_sample_std)
p <- dim(z_sample_std)[2]


### 2. Initialize objects to save estimation results
iter_n <- length( num_knot )

results_list <- vector(mode = "list", length = iter_n)

est_K0 <- rep(0, iter_n)
est_eta <- matrix(0, iter_n, p)
est_beta <- matrix(0, iter_n, p)
est_lambda0 <- matrix(0, iter_n, p)
est_lambda1 <- matrix(0, iter_n, p)
est_lambda2 <- matrix(0, iter_n, p)




set.seed(12345)

### 3. Start for loop for each number of knots
for(i in 1:iter_n){

  cat(i, "th iteration begun \n")

  ### 4. Compute B_tilde matrix
  knots_set <- (num_knot - 7)[i]
  knots_vec <- select.knots(t_vec, knots = knots_set ,p=3)
  B_mat <- splineDesign(knots_vec, t_vec, ord = 4)

  G <- (t(B_mat) %*% B_mat / m)
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
  fit_face <- fpca_face_adjust(Y=y_demeaned_mat,ydata=NULL,Y.pred = NULL,argvals=t_vec,pve = 0.99, npc  = NULL,
                             var = TRUE, simul = FALSE, sim.alpha = 0.95,
                             center=FALSE,knots_para=knots_set,p=3,m=2,lambda=NULL,alpha = 1,
                             search.grid=TRUE,search.length=100,
                             method="L-BFGS-B", lower=-20,upper=20, control=NULL)
  sigma2_hat <- fit_face$sigma2
  cov_y_face <- fit_face$Cov
  cov_y_sample <-  1/(m) * t(y_demeaned_mat) %*% (y_demeaned_mat) / (n-1)


  ### 7. CoFA for common component part
  CoFA_fit <- CoFA_CrossCov(Y_cent = y_demeaned, Z_cent = z_sample_std_mat, B_tilde = B_tilde, gamma = 2,
                            control = list(num_rep=100, num_fold=5, num_tau=100))
  cov_y <- t(y_demeaned_mat) %*% (y_demeaned_mat) / (n-1)
  C_tilde <- 1/(m^2) * t(B_tilde) %*% ( cov_y - diag(rep(sigma2_hat,m)  ) ) %*% B_tilde

  K0_est <- CoFA_fit$rank
  est_K0[i] <- K0_est
  if(K0_est == 1){
    est_eta[i, 1] <- CoFA_fit$D
    est_beta[i,1] <- est_eta[i, 1] / (t(CoFA_fit$U) %*% C_tilde %*% CoFA_fit$U   )
    est_lambda0[i, 1] <-  (t(CoFA_fit$U) %*% C_tilde %*% CoFA_fit$U  )

  } else {

    est_eta[i,1:length(diag(CoFA_fit$D))] <- diag(CoFA_fit$D)
    for(j in 1:K0_est){
      est_beta[i,j] <- diag(CoFA_fit$D)[j] / (t(CoFA_fit$U[,j]) %*% C_tilde %*% CoFA_fit$U[,j]   )
      est_lambda0[i, j] <-  (t(CoFA_fit$U[,j]) %*%  C_tilde %*% CoFA_fit$U[,j]  )
    }
  }




  ### 8. Estimation for independent component part
  ## 8-a functional data
  if(K0_est == 1){

    A0 <- CoFA_fit$U
    P1 <- diag(1, dim(B_mat)[2]) - A0 %*% solve( t(A0) %*% A0) %*% t(A0)
    svd_k1 <- svd( P1 %*% C_tilde  %*% t(P1) )

  } else{

    A0 <- CoFA_fit$U[,1:K0_est]
    P1 <- diag(1, dim(B_mat)[2]) - A0 %*% solve( t(A0) %*% A0) %*% t(A0)
    svd_k1 <- svd( P1 %*% C_tilde  %*% t(P1) )

  }
  est_lambda1[i,] <- svd_k1$d[1:p]


  ## 8-b multivariate data
  Cov_Z <- cov(z_sample_std_mat)
  V0 <- CoFA_fit$V
  P_v0 <- V0  %*% solve( t(V0) %*% V0  )%*% t(V0)

  # estimation for sigma^{2}_{omega}
  # based on Ordinary Least Square with a constraint
  var_omega <- var_omega_est(z_data_cent = z_sample_std_mat, V0_mat = CoFA_fit$V,
                             est_beta_vec =  est_beta[i, 1:K0_est],
                             est_eta_vec = est_eta[i, 1:K0_est],
                             est_K0 = K0_est)
  error_omega_est <-  var_omega
  # end of sigma^{2}_{omega} estimation

  P_v1 <- diag(1, p) - P_v0
  svd_nu_k1 <- svd(P_v1 %*% (Cov_Z - diag(error_omega_est , p)) %*% P_v1)
  est_lambda2[i, ] <- svd_nu_k1$d



  ### 9. mBLUP
  ## 9-a mMLUP for common components
  new_data <- rbind( 1/sqrt(m) * t(y_sample), t(z_sample_std_mat) )
  mean_data <- rbind(1/sqrt(m) * t(y_centered), matrix(0, p, n) )
  x_demeaned <- ( new_data - mean_data  )

  cov_x_demeaned <- x_demeaned %*% t(x_demeaned) / (n-1)
  if(K0_est ==1){
    M_3 <- sqrt( 1/(1+(est_beta[i,1])^2 ) ) * cbind( 1/sqrt(m) * t(B_tilde %*% CoFA_fit$U ), (est_beta[i,1]) * t(CoFA_fit$V) )
    x_new <- M_3 %*% x_demeaned
    Blup_pred_com <-  est_lambda0[i,1] *  sqrt(1 + (est_beta[i, 1:K0_est])^2) * solve( M_3 %*% cov_x_demeaned %*% t(M_3) ) %*% x_new
  } else {
    M_3 <- diag( sqrt( 1/(1+(est_beta[i,1:K0_est])^2 )  ) ) %*% cbind(  1/sqrt(m) * t( B_tilde %*% CoFA_fit$U[, 1:K0_est] ),
                                                                        diag( (est_beta[i,1:K0_est]) ) %*% t(CoFA_fit$V[,1:K0_est]))
    x_new <- M_3 %*% x_demeaned
    Blup_pred_com <- diag( est_lambda0[i,1:K0_est]) %*% diag( sqrt(1 + (est_beta[i, 1:K0_est])^2) )   %*% solve( M_3 %*% cov_x_demeaned %*% t(M_3) ) %*% x_new
  }


  ## 9-b mMLUP for independent components in functional data
  y_demeaned_blup <- 1/sqrt(m) * ( t(y_sample) - t(y_centered)  )
  Phi_1 <- 1/sqrt(m) * B_tilde %*% svd_k1$u[,1:3]
  Blup_pred_ind_func <- diag( svd_k1$d[1:3] )   %*% solve( t(Phi_1) %*% y_demeaned_blup %*% t(y_demeaned_blup) %*% Phi_1 / (n-1) ) %*% t(Phi_1) %*% y_demeaned_blup


  ## 9-c mMLUP for independent components in multivariate data
  Phi_2 <- svd_nu_k1$u[,1:3]
  Blup_pred_ind_multi <- diag( (svd_nu_k1$d[1:3] ) )   %*% solve(t(Phi_2) %*% Cov_Z %*% Phi_2 ) %*% t(Phi_2) %*% t(z_sample_std_mat)


  ### 10. Save all results
  results_list[[i]] <- list( "CoFA_fit"=CoFA_fit,
                             "svd_k1"=svd_k1, "svd_nu_k1"=svd_nu_k1,
                             "cov_y_face"=cov_y_face, "Cov_Z"=Cov_Z,
                             "fit_face"=fit_face,
                             "cov_y_sample"=cov_y_sample,
                             "y_centered"= y_centered,
                             "Blup_pred_com"=Blup_pred_com,
                             "Blup_pred_ind_func"=Blup_pred_ind_func,
                             "Blup_pred_ind_multi"=Blup_pred_ind_multi,
                             "B_mat"=B_mat, "B_tilde"= B_tilde)
}

if(which_ear == "right"){
  # save(results_list, file="./Data/real_data_app_right.RData")
  save(results_list, z_sample_std_mat, file="real_data_app_right.RData")

} else if(which_ear == "left"){
  # save(results_list, file="./Data/real_data_app_left.RData")
  save(results_list, z_sample_std_mat, file="real_data_app_left.RData")

}

res <- list("results_list"=results_list, "z_sample_std_mat"=z_sample_std_mat)

return(res)
}


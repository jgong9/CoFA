#' Common Factor Analysis between Functional and Multivariate Data by Adaptive Nuclear Norm Penalty
#'
#' @export CoFA_CrossCov
#' @param Y_cent The centered functional data matrix
#' @param Z_cent The centerd multivariate data matrix
#' @param t_vec A vector of the observation points of functional data
#' @param B_tilde The orthogonal B-spline basis matrix
#' @param gamma The parameter used for weights in adaptive nuclear norm penalty
#' @param control A list that contains the cross validation setting
#' @return A list that contains the decomposition of estimated M matrix with its rank and the chosen tuning parameter

CoFA_CrossCov <- function(Y_cent , Z_cent, t_vec, B_tilde, gamma=2,
                          control = list(num_rep=100, num_fold=5, num_tau=100)
                          ){

                            m <- dim(Y_cent)[2]
                            n_original <- nrow(Y_cent)

                            Cross_cov <- 1/(n_original-1) *( t(Y_cent) %*% (Z_cent) )

                            w <- quadWeights(argvals = t_vec , method = "trapezoidal")
                            W_mat <- diag(w)
                            Gamma_mat <- t(B_tilde) %*% W_mat %*% Cross_cov


                            svd_Gamma <- svd(Gamma_mat)
                            d_origin <- svd_Gamma$d
                            d_result <- vector("numeric", length(d_origin))
                            est_rank <- 0
                            max_rank <- min(dim(Gamma_mat)[1], dim(Gamma_mat)[2])

                            Adx <- (svd_Gamma$d) ^ (1+gamma)

                            tau_vec <- exp(seq(log(max(Adx)), log(Adx[max_rank] ),
                                                  length = control$num_tau) )

                            #######################################################################
                            #### Cross Validation

                              tau_final_vec <- rep(0, control$num_rep)

                              #### (1) Repeated CV for num_rep times without seed number
                              for(loc_i in 1:control$num_rep){

                                order_new <- sample(1:n_original, n_original,replace=F)
                                if(n_original %%control$num_fold != 0 ){
                                  number_per_fold <- n_original %/%control$num_fold
                                  add_last_fold <- n_original %%control$num_fold

                                  label_vec <- rep(0, n_original)
                                  for( k in 1:control$num_fold){
                                    label_vec[order_new[((k-1)*(number_per_fold) +1 ):( k * number_per_fold ) ]] <- k
                                  }

                                  label_vec[  order_new[  (control$num_fold * number_per_fold + 1):(control$num_fold * number_per_fold + add_last_fold  ) ] ] <-control$num_fold
                                } else {
                                  number_per_fold <- n_original %/%control$num_fold

                                  label_vec <- rep(0, n_original)
                                  for( k in 1:control$num_fold){
                                    label_vec[order_new[((k-1)*(number_per_fold) +1 ):( k * number_per_fold ) ]] <- k
                                  }
                                }


                                sse_CV_mat <- matrix(0, control$num_tau,control$num_fold)

                                ### (2) Start of for loop at "loc_i"th iteration

                                for(l in 1:control$num_tau){
                                  tau_tun <- tau_vec[l]


                                  ## (3) Start of for loop with a given tau
                                  for(k_cv in 1:control$num_fold){
                                    Y_cent_train <- Y_cent[ label_vec != k_cv, ]
                                    Y_cent_test <- Y_cent[label_vec == k_cv , ]

                                    Z_cent_train <- Z_cent[ label_vec != k_cv, ]
                                    Z_cent_test <- Z_cent[label_vec == k_cv , ]



                                    Cross_cov_train <- 1/(nrow(Y_cent_train) - 1) * (t(Y_cent_train) %*% Z_cent_train)
                                    Gamma_mat_train <- t(B_tilde) %*% W_mat %*% Cross_cov_train

                                    Cross_cov_test <- 1/(nrow(Y_cent_test) - 1) * (t(Y_cent_test) %*% Z_cent_test)
                                    Gamma_mat_test <- t(B_tilde) %*% W_mat %*% Cross_cov_test


                                    svd_Gamma_train <- svd(Gamma_mat_train)
                                    d_origin <- svd_Gamma_train$d
                                    d_result <- vector("numeric", length(d_origin))
                                    est_rank <- 0

                                    shrinked <- d_origin - tau_tun * ( d_origin^(- gamma ) )

                                    if(sum(shrinked < 0) > 0){
                                      est_rank <- which(shrinked < 0 )[1] - 1
                                    } else {
                                      est_rank <- length(d_origin)
                                    }


                                    if(est_rank > 0 ){
                                      ## D for different rank
                                      if(est_rank == 1){
                                        D = shrinked[1:est_rank]
                                      } else {
                                        D = diag( shrinked[1:est_rank] )
                                      }
                                      ##

                                      res <- list(
                                        U = svd_Gamma_train$u[,1:est_rank],
                                        V = svd_Gamma_train$v[,1:est_rank],
                                        D= D,
                                        rank=est_rank
                                      )
                                    } else {
                                      res <- list(rank=est_rank)
                                    }


                                    if(res$rank != 0){
                                      if(res$rank ==1 ){
                                        temp_M <-  res$D * res$U %*% t(res$V)
                                      } else {
                                        temp_M <- res$U %*% res$D %*% t(res$V)
                                      }
                                    } else {
                                      ## Null solution
                                      temp_M <- matrix(0, dim(Gamma_mat_test)[1], dim(Gamma_mat_test)[2])
                                    }
                                    sse_CV_mat[l, k_cv] <- sum( (Gamma_mat_test - temp_M)^2 )


                                  } ## (3) End of for loop with a given tau
                                } ### (2) End of for loop for all taus

                                # Compute Mean of 5-fold test error for each lambda candidate
                                # and choose the one which achieved the minimum average errors.
                                CV_tau <- tau_vec[which.min(rowMeans(sse_CV_mat))]

                                tau_final_vec[loc_i] <- CV_tau
                              } # (1) End of all repetitions

                              ###### 4. Take the average of tau_final_vec as the final tuning value
                              tau_final <- mean(tau_final_vec)


                              ### End of CV process    ###################################################
                              ############################################################################

                              svd_Gamma <- svd(Gamma_mat)
                              d_origin <- svd_Gamma$d
                              est_rank <- 0

                              shrinked <- d_origin - tau_final * ( d_origin^(-2) )


                              if(sum(shrinked < 0) > 0){
                                est_rank <- which(shrinked < 0 )[1] - 1
                              } else {
                                est_rank <- length(d_origin)
                              }

                              if(est_rank > 0 ){
                                ## D for different rank
                                if(est_rank == 1){
                                  D = shrinked[1:est_rank]
                                } else {
                                  D = diag( shrinked[1:est_rank] )
                                }
                                ##

                                res <- list(
                                  U = as.matrix(svd_Gamma$u[,1:est_rank]),
                                  V = as.matrix(svd_Gamma$v[,1:est_rank]),
                                  D = D,
                                  rank=est_rank,
                                  tau_final = tau_final,
                                  tau_final_vec = tau_final_vec
                                )
                              } else {
                                res <- list(
                                            U = NULL,
                                            V = NULL,
                                            D = NULL,
                                            rank=est_rank,
                                            tau_final = tau_final,
                                            tau_final_vec = tau_final_vec
                                            )
                              }

                              return(res)


                          }



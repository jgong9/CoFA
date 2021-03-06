% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CoFA_CrossCov.R
\name{CoFA_CrossCov}
\alias{CoFA_CrossCov}
\title{Common Factor Analysis between Functional and Multivariate Data by Adaptive Nuclear Norm Penalty}
\usage{
CoFA_CrossCov(
  Y_cent,
  Z_cent,
  t_vec,
  B_tilde,
  gamma = 2,
  control = list(num_rep = 100, num_fold = 5, num_tau = 100)
)
}
\arguments{
\item{Y_cent}{The centered functional data matrix}

\item{Z_cent}{The centerd multivariate data matrix}

\item{t_vec}{A vector of the observation points of functional data}

\item{B_tilde}{The orthogonal B-spline basis matrix}

\item{gamma}{The parameter used for weights in adaptive nuclear norm penalty}

\item{control}{A list that contains the cross validation setting}
}
\value{
A list that contains the decomposition of estimated M matrix with its rank and the chosen tuning parameter
}
\description{
Common Factor Analysis between Functional and Multivariate Data by Adaptive Nuclear Norm Penalty
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_study.R
\name{simulation_study}
\alias{simulation_study}
\title{A function to conduct simulation study at a given signal noise level and a number of iterations
Note that two seed numbers are implanted in this function: the one for the true eigen vectors, and the other for reproducibility of data generation}
\usage{
simulation_study(
  snr = 5,
  sample_size = c(100, 200, 500),
  N = 1000,
  num_resampling_CV = 100
)
}
\arguments{
\item{snr}{The signal noise ratio value}

\item{sample_size}{A numeric vector of sample sizes}

\item{N}{The number of iterations}

\item{num_resampling_CV}{The number of repetitions for CoFA function}
}
\value{
A list that contains 3 lists as elements for 3 sample sizes (100, 200, 500): each of them has N+1 elements for the first N having CoFA fits and BLUPS at a single iteration and the last one including the results of parameter estimation and assessment
}
\description{
A function to conduct simulation study at a given signal noise level and a number of iterations
Note that two seed numbers are implanted in this function: the one for the true eigen vectors, and the other for reproducibility of data generation
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_real_data.R
\name{read_real_data}
\alias{read_real_data}
\title{A function to read in the real data of NHANES}
\usage{
read_real_data(select_ear = "right")
}
\arguments{
\item{select_ear}{Either "right" or "left" to choose which ear of data is loaded}
}
\value{
A list that contains two data: the one for functional and multivariate and the other for the label of tympnometry
}
\description{
A function to read in the real data of NHANES
}

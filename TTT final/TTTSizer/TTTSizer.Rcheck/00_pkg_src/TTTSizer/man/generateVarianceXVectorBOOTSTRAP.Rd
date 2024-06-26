% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools.R
\name{generateVarianceXVectorBOOTSTRAP}
\alias{generateVarianceXVectorBOOTSTRAP}
\title{TTTSizer}
\usage{
generateVarianceXVectorBOOTSTRAP(x, bootFactor = 100)
}
\arguments{
\item{x}{The original vector with your data points.}

\item{bootFactor}{The number of bootstrap iterations, default is 100}
}
\value{
A matrix with all the coovariances and variances. The main diagonal contains the variance vector
}
\description{
Generate the Variance/Coovariance matrix with a bootstrap algorithm
}

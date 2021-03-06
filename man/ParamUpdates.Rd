% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{ParamUpdates}
\alias{ParamUpdates}
\title{This function updates the block-specific parameters.}
\usage{
ParamUpdates(betar, etar, xr, yr, dat, ur, rr, beta, rho, alpha, tau, n, ni)
}
\arguments{
\item{betar}{a vector containing the current values for the block-specific beta parameters}

\item{etar}{a vector containing the current values for the block-specific eta parameters}

\item{xr}{a matrix containing the block-specific design matrix}

\item{yr}{a matrix containing the block-specific outcome vector}

\item{dat}{a matrix containing the block-specific matrix inverse of t(xr)%*%xr + I}

\item{ur}{a vector containing the current values for the block-specific u variables}

\item{beta}{a vector containing the current values for the global beta parameters}

\item{rho}{a numeric scalar}

\item{alpha}{a numeric scalar}

\item{tau}{a numeric scalar}

\item{n}{an integer giving the overall sample size for the full data}

\item{ni}{an integer giving the block-specific sample size}
}
\description{
This function updates the block-specific parameters.
}

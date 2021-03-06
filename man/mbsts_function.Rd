% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mbsts_function.R
\docType{methods}
\name{mbsts_function}
\alias{mbsts_function}
\alias{mbsts_function,array-method}
\title{Main function for the multivariate Bayesian structural time series (MBSTS) model}
\usage{
mbsts_function(
  Y,
  Xtrain,
  STmodel,
  ki,
  pii,
  b = NULL,
  v0,
  kapp = 0.01,
  R2 = 0.8,
  v = 0.01,
  ss = 0.01,
  mc = 500,
  burn = 50
)

\S4method{mbsts_function}{array}(
  Y,
  Xtrain,
  STmodel,
  ki,
  pii,
  b = NULL,
  v0,
  kapp = 0.01,
  R2 = 0.8,
  v = 0.01,
  ss = 0.01,
  mc = 500,
  burn = 50
)
}
\arguments{
\item{Y}{A (\eqn{n*m})-dimensional matrix containing multiple target series, where \eqn{n} is the number of observations and \eqn{m} is the number of target series.}

\item{Xtrain}{A (\eqn{n*K})-dimensional matrix containing all candidate predictor series for each target series. \eqn{K=\sum k_i} is the number of all candidate predictors for all target series. The first \eqn{k_1} variables are the set of candidate predictors for the first target series, and the next \eqn{k_2} variables are the set of candidate predictors for the second target series, etc. Note that, one variable can appear in the X.star several times, since different target series can contain the same candidate predictors.}

\item{STmodel}{A state space model of SSmodel class returned by tsc.setting.}

\item{ki}{A vector of integer values denoting the acumulated number of predictors for target series. For example, if there are three target series where the first has \eqn{8} predictors, the second has \eqn{6} predictors, and the third has \eqn{10} predictors, then the vector is c(\eqn{8,14,24}).}

\item{pii}{A vector describing the prior inclusion probability of each candidate predictor.}

\item{b}{NULL or a vector describing the prior means of regression coefficients. The default value is NULL.}

\item{v0}{A numerical value describing the prior degree of freedom of the inverse Wishart distribution for \eqn{\Sigma_\epsilon}.}

\item{kapp}{A scalar value describing the number of observations worth of weight on the prior mean vector. The default value is \eqn{0.01}.}

\item{R2}{A numerical value taking value in \eqn{[0,1]}, describing the expected percentage of variation of \eqn{Y} to be explained by the model. The default value is \eqn{0.8}.}

\item{v}{A numerical value describing the prior degree of freedom of the inverse Wishart distribution for (\eqn{\Sigma_\mu,\Sigma_\delta,\Sigma_\tau,\Sigma_\omega}). The default value is \eqn{0.01}.}

\item{ss}{A numerical value describing the prior scale matrix of the inverse Wishart distribution for (\eqn{\Sigma_\mu,\Sigma_\delta,\Sigma_\tau,\Sigma_\omega}). The default value is \eqn{0.01}.}

\item{mc}{A positive integer giving the desired number of MCMC draws. The default value is \eqn{500}.}

\item{burn}{A positive integer giving the number of initial MCMC draws to be discarded. The default value is \eqn{50}.}
}
\value{
An object of mbsts class
}
\description{
The MBSTS model uses MCMC to sample from the posterior distribution of a MBSTS model. The model is given by
\deqn{y=\mu+\tau+\omega+\beta X+\epsilon,}
where \eqn{\mu}, \eqn{\tau}, \eqn{\omega}, \eqn{\beta X}, and \eqn{\epsilon} denote the trend component, the seasonal component, the cycle component, the regression component, and the error term, respectively. Note that, without a regression component, the MBSTS model is an ordinary state space time series model. The predictors and response variables in the MBSTS model are designed to be contemporaneous. Lags and differences can be generated by manipulating the predictor matrix. The "spike-and-slab" prior is used for the regression component of models, which enables feature selection among a large number of features.
}
\references{
\Qiu2018

\Ning2021

\Jammalamadaka2019
}
\author{
Jinwen Qiu \email{qjwsnow_ctw@hotmail.com} Ning Ning \email{patricianing@gmail.com}
}

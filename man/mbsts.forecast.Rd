% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mbsts.forecast.R
\name{mbsts.forecast}
\alias{mbsts.forecast}
\title{Specification of time series components}
\usage{
mbsts.forecast(object, STmodel, newdata, steps = 1)
}
\arguments{
\item{object}{An object of the mbsts class created by a call to the mbsts_function function.}

\item{STmodel}{An object of the SSModel class created by a call to the tsc.setting function.}

\item{newdata}{A vector or matrix containing the predictor variables to use in making the prediction. This is only required if the mbsts model has a regression component.}

\item{steps}{An integer value describing the number of time steps ahead to be forecasted. If it is greater than the number of new observations in the newdata, zero values will fill in missing new observations.}
}
\value{
An object of predicted values which is a list containing the following:
\item{pred.dist}{An array of draws from the posterior predictive distribution. The first dimension in the array represents time, the second dimension denotes each target series, and the third dimension indicates each MCMC draw.}
\item{pred.mean}{A matrix giving the posterior mean of the prediction for each target series.}
\item{pred.sd}{A matrix giving the posterior standard deviation of the prediction for each target series.}
\item{pred.se}{A matrix giving the posterior standard error of the prediction for each target series, calculated by pred.sd divided by the square root of the numer of MCMC iterations.}
}
\description{
Generate draws from the posterior predictive distribution of a mbsts object. Samples from the posterior predictive distribution of the MBSTS model.
}
\references{
\Qiu2018

\Ning2021

\Jammalamadaka2019
}
\author{
Jinwen Qiu \email{qjwsnow_ctw@hotmail.com} Ning Ning \email{patricianing@gmail.com}
}

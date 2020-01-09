\name{tsc.setting}
\alias{tsc.setting}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Specification of time series components
}
\description{
Specify three time series components for MBSTS model (generalized linear trend, seasonality and cycle). A generalization of the local linear trend model where the slope exhibits stationarity instead of obeying a random walk, is expressed in the form as:
\deqn{\mu_{t+1}=\mu_t+\delta_t+rnorm(0,\Sigma_\mu)}
\deqn{\delta_{t+1}=D+\rho(\delta_t-D)+rnorm(0,\Sigma_\delta)}
The seasonal model in the time domain is:
\deqn{\tau_{t+1}^{(i)}=-\sum_{k=0}^{S_i-2}\tau ^{(i)}_{t-k}+rnorm(0,\Sigma_\tau)}
The cycle component is postulated as:
\deqn{\omega_{t+1}=varrho*cos(\lambda)\omega_{t}+varrho*sin(\lambda)\omega*_{t}+rnorm(0,\Sigma_\omega)}
\deqn{\omega*_{t+1}=-varrho*sin(\lambda)\omega_{t}+varrho*cos(\lambda)\omega*_{t}+rnorm(0,\Sigma_\omega)}
}
\usage{
tsc.setting(Y, mu = NULL, rho = NULL, S = NULL, vrho = NULL, lambda = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
The multivariate time series to be modeled, as a numeric matrix convertible to xts. This state model assumes Y contain daily data.
}
  \item{mu}{
A vector of logic values, indicating whether to include a local trend for each target series.
}
  \item{rho}{
A vector of numerical values, all of which are between 0 and 1, describing the learning rates at which the local trend is updated for each target series. The value 0 in jth entry indicates that jth target series does not include slope of trend. 
}
  \item{S}{
A vector of integer values, representing the number of seasons to be modeled for each target series. The value 0 in jth entry indicates that jth target series does not include seasonal component. 
}
  \item{vrho}{
A vector of numerical values between 0 and 1, describing a damping factor for each target series. The value 0 in jth entry indicates that jth target series does not include cycle component. 
}
  \item{lambda}{
A vector of numerical values, each of them is equal to \eqn{2\pi/q}, describing 
the frequency with q being a period such that \eqn{0<\lambda<\pi}.
}
}

\details{
The function allows users to include different combination of time series components for each target series. 
}

\value{
Returns a customerized State Space model, which is an object of class SSModel.
}

\references{
Harvey (1990), "Forecasting, structural time series, and the Kalman filter", Cambridge University Press.

Durbin and Koopman (2001), "Time series analysis by state space methods", Oxford University Press.
}
\author{
Jinwen Qiu <jqiu@pstat.ucsb.edu>
}
\seealso{
\code{\link[KFAS]{SSModel}}
}
\examples{
data(exdata)
Y<-as.matrix(exdata[,1:2])
STmodel<-tsc.setting(Y,mu=c(1,1),rho=c(0.6,1),S=c(4,0),
                     vrho=c(0,0.5),lambda=c(0,pi/10))
}
\keyword{state components}% use one of  RShowDoc("KEYWORDS")


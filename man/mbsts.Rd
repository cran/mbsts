\name{mbsts}
\alias{mbsts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate Bayesian Structural Time Series
}
\description{
Uses MCMC to sample from the posterior distribution of a Multivariate Bayesian structural time series model. This function can be used either with or without contemporaneous predictor variables (in a time series regression).

If predictor variables are present, the regression coefficients are fixed. The predictors and response in the formula are contemporaneous, so if you want lags and differences you need to put them in the predictor matrix yourself.
If no predictor variables are used, then the model is an ordinary state space time series model.

The model allows for several useful extensions beyond standard Bayesian dynamic linear models. A spike-and-slab prior is used for the (static) regression component of models that include predictor variables. This is especially useful with large numbers of regressor series. 

The model in state space form can be written as:
\deqn{y_t=\mu_t+\tau_t+\omega_t+\beta*X_t+rnorm(0,\Sigma_\epsilon)}
where \eqn{\mu},\eqn{\tau} and \eqn{\omega} denotes trend, seasonal and cycle component, respectively.
}
\usage{
mbsts(Y, 
      X.star = NULL, 
      STmodel = NULL, 
      ki = NULL, 
      pii = NULL, 
      b = NULL, 
      kapp = 0.1, 
      R2 = 0.8, 
      v0 = NULL, 
      v = 0.01, 
      ss = 0.01, 
      mc = 500, 
      burn = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
An n by m matrix contains multiple target series. The n and m is the number of observations and target series, respectively.
}
  \item{X.star}{
An n by K matrix contains all candidate predictor series for each target series. 
The K=\eqn{\sum k_i} is the number of all candidate predictors for all target series. The first \eqn{k_1} variables make up the set of candidate predictors for 1st target series, and the next \eqn{k_2} variables are candidate predictors for 2nd target series, etc.
One variable can appear in the X.star several times , since different target series can contain same candidate predictors.
}
  \item{STmodel}{
A state space model of SSmodel object returned by tsc.setting. 
}
  \item{ki}{
A vector of integer values, denoting the location of last candidate predictor for each target series, such as \eqn{c(k_1,k_1+k_2,...,K)}
}
  \item{pii}{
A column matrix with length equal to K, describing the prior inclusion probability of each candidate predictor. 
}
  \item{b}{
A column matrix, describing the prior means for regression coefficients. 
}
  \item{kapp}{
A scalar value, describing the number of observations worth of weight on the prior mean vector.
}
  \item{R2}{
A numerical value from 0 to 1, describing the expected percent of variation of Y to be explained by model.
}
  \item{v0}{
A numerical value, describing the degree of freedom for prior Inverse Wishart distribution (\eqn{\Sigma_\epsilon}).
}
  \item{v}{
A numerical value, describing the prior sample size of Inverse Gama distribution (\eqn{\Sigma_\mu,\Sigma_\delta,\Sigma_\tau,\Sigma_\omega}).
}
  \item{ss}{
A numerical value, describing the prior sum of squares of Inverse Gama distribution (\eqn{\Sigma_\mu,\Sigma_\delta,\Sigma_\tau,\Sigma_\omega}).
}
  \item{mc}{
A positive integer giving the desired number of MCMC draws.
}
  \item{burn}{
A positive integer giving the initial number of MCMC draws to be discarded.
}
}
\details{
See the paper worte by Jammalamadaka, Qiu and Ning (2018) <https://arxiv.org/pdf/1801.03222.pdf> 
}
\value{
An object of class mbsts which is a list with the following components:
  \item{Ind}{One K by mc-burn matrix, containing MCMC draws of indicator variable (gama). If the X.star is null, it will not be returned.
}
  \item{beta.hat}{
One K by mc-burn matrix, containing MCMC draws of regression coefficients (beta).If the X.star is null, it will not be returned.
}
  \item{B.hat}{This is beta.hat in matrix form. If the X.star is null, it will not be returned.
}
  \item{ob.sig2}{
An array with dimension (m*m*(mc-burn)),containing MCMC draws of variance and covariance matrix for residuals. 
}
  \item{States}{
An array with dimension (n*m1*(mc-burn)),containing MCMC draws of all time series components. The  m1 is the number of all time series components, respectively. If the STmodel is null, it will not be returned.
}
  \item{st.sig2}{
One K by mc-burn matrix, containing MCMC draws of variance for time series components. If the STmodel is null, it will not be returned.
}
}
\references{
Harvey (1990),"Forecasting, structural time series, and the Kalman filter", Cambridge University Press.

Durbin and Koopman (2001), "Time series analysis by state space methods", Oxford University Press.

Goerge and McCulloch (1997) "Approaches for Bayesian variable selection", Statistica Sinica pp 339-74.

Scott and Varian (2014) "Predicting the present with Bayesian structural time series", International Journal of Mathematical Modelling and Numerical Optimisation, 5 (1-2), 4-23.
  
Griffiths (2003) "Bayesian inference in the seemingly unrelated regressions model",
Computer aided econometrics, pp 263-290, CRC Press.
}
\author{
Jinwen Qiu <jqiu@pstat.ucsb.edu>
}
\seealso{
\code{\link{tsc.setting}}
}
\examples{
data(exdata)

#Two target series
Y<-as.matrix(exdata[,1:2])
#Sixteen candidate predictors
X.star<-as.matrix(exdata[,3:18])

#split dataset into training set and test set
n=dim(Y)[1]
ntrain=n-5
Ytrain<-Y[1:ntrain,]
Xtrain<-X.star[1:ntrain,]
Ytest<-Y[(ntrain+1):n,]
Xtest<-X.star[(ntrain+1):n,]

#Specify time series components
STmodel<-tsc.setting(Ytrain,mu=c(1,1),rho=c(0.6,1),S=c(4,0),
                     vrho=c(0,0.5),lambda=c(0,pi/10))
                     
#prior parameters setting 
#gama
ki<- c(8,dim(Xtrain)[2])
pii<- matrix(rep(0.5,dim(Xtrain)[2]),nrow=dim(Xtrain)[2])

#beta
b<-matrix(0,dim(Xtrain)[2])
kapp<-0.01

#v0 and V0 for obs Sigma
R2<-0.8
v0<-5

#State component Sigma
v<-0.01
ss<-0.01

#train a mbsts model
mbsts.model<-mbsts(Ytrain,Xtrain,STmodel,ki,pii,b,kapp,R2,v0,v,ss,mc=15,burn=5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{model}% use one of  RShowDoc("KEYWORDS")


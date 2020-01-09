\name{mbsts.forecast}
\alias{mbsts.forecast}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multi-steps ahead Forecast by Multivariate Bayesian Structural Time Series Model
}
\description{
Generated draws from the posterior predictive distribution of a mbsts object.
}
\usage{
mbsts.forecast(mbsts, STmodel = NULL, newdata = NULL, steps = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mbsts}{
An object of class mbsts created by a call to the function mbsts.
}
  \item{STmodel}{
An object of class SSModel created by a call to the function tsc.setting.
}
  \item{newdata}{
a vector or matrix containing the predictor variables to use in making the prediction. This is only required if object contains a regression component.
}
  \item{steps}{
An integer value, describing the number of time steps ahead to be forecasted. If it is greater than the number of new observations in the newdata, the zero values will fill in missing new observations.
}
}
\details{
Samples from the posterior distribution of a Bayesian structural time series model. This function can be used either with or without contemporaneous predictor variables (in a time series regression).

If predictor variables are present, the regression coefficients are fixed
The predictors and response in the formula are contemporaneous, so if you want lags and differences you need to put them in the predictor matrix yourself.

If no predictor variables are used, then the model is an ordinary state space time series model.
}
\value{
An object of predicted values by mbsts which is a list with the following components:
\item{pred.distribution}{An array of draws from the posterior predictive distribution. The first dimension in the array represents time, the second dimension denotes each target series, and the third dimension indicates each MCMC draw.}
\item{pred.mean}{A matrix giving the posterior mean of the prediction for each target series.}

}

\references{
Harvey (1990),"Forecasting, structural time series, and the Kalman filter", Cambridge University Press.

Durbin and Koopman (2001), "Time series analysis by state space methods", Oxford University Press.

Goerge and McCulloch (1997) "Approaches for Bayesian variable selection", Statistica Sinica pp 339-74.

Scott and Varian (2014) "Predicting the present with Bayesian structural time series", International Journal of Mathematical Modelling and Numerical Optimisation, 5 (1-2), 4-23.
}
\author{
Jinwen Qiu <jqiu@pstat.ucsb.edu>
}
\seealso{
\code{\link{mbsts}}
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

#make a 5-steps prediction
output<-mbsts.forecast(mbsts.model,STmodel,newdata=Xtest,steps=3)
#error<-abs(output$pred.mean-Ytest)
#error
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{forecast}% use one of  RShowDoc("KEYWORDS")


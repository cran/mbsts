## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(mbsts)

## -----------------------------------------------------------------------------
###############Setup###########
n<-505 #n: sample size
m<-2 #m: dimension of target series

cov<-matrix(c(1.1,0.7,0.7,0.9), nrow=2, ncol=2) #covariance matrix of target series 

###############Regression component###########
#coefficients for predictors
beta<-t(matrix(c(2,-1.5,0,4,2.5,0,0,2.5,1.5,-1,-2,0,0,-3,3.5,0.5),nrow=2,ncol=8)) 

set.seed(100)
#predictors
X1<-rnorm(n,5,5^2)
X4<-rnorm(n,-2,5)
X5<-rnorm(n,-5,5^2)
X8<-rnorm(n,0,100)
X2<-rpois(n, 10)
X6<-rpois(n, 15)
X7<-rpois(n, 20)
X3<-rpois(n, 5)
X<-cbind(X1,X2,X3,X4,X5,X6,X7,X8) 

###############Simulated data################
set.seed(100)
data=sim_data(X=X, beta=beta, cov, k=c(8,8), mu=c(1,1), rho=c(0.6,0.8), 
              Dtilde=c(-1,3), Season=c(100,0),      
              vrho=c(0,0.99),lambda=c(0,pi/100))


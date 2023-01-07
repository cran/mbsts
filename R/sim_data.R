#' @title Simulate data
#' @author Jinwen Qiu \email{qjwsnow_ctw@hotmail.com} Ning Ning \email{patricianing@gmail.com}
#' @description Generate simulated data in the form of structural time series
#' @importFrom stats rnorm
#' @docType methods

#' @param X  A (\eqn{n*K})-dimensional matrix containing predictors, where \eqn{n} is the number of observations. \eqn{K=\sum k_i} is the number of all candidate predictors for all target series. The first \eqn{k_1} variables are the set of candidate predictors for the first target series, and the next \eqn{k_2} variables are the set of candidate predictors for the second target series, etc. 
#' @param beta  A (\eqn{K*m})-dimensional matrix containing all candidate predictor series for each target series. 
#' @param cov  A (\eqn{m*m})-dimensional matrix containing covariances
#' @param k  A \eqn{m}-dimensional array containing the number of candidate predictors for each of the \eqn{m} target series.
#' @param mu A \eqn{m}-dimensional array with \eqn{1} representing modeling with trend for this target time series.
#' @param rho  A \eqn{m}-dimensional array representing the learning rates at which the local trend is updated.
#' @param mean_trend  A numerical value standing for the mean of the error term of the trend component. The default value is \eqn{1}.
#' @param sd_trend  A numerical value standing for the standard deviation of the error term of the trend component. The default value is \eqn{0.5}.
#' @param mean_season A numerical value standing for the mean of the error term of the seasonal component. The default value is \eqn{20}.
#' @param sd_season A numerical value standing for the standard deviation of the error term of the seasonal component. The default value is \eqn{0.5}.
#' @param mean_cycle  A numerical value standing for the mean of the error term of the cycle component. The default value is \eqn{20}.
#' @param sd_cycle  A numerical value standing for the standard deviation of the error term of the cycle component. The default value is \eqn{0.5}.
#' @param Dtilde  A \eqn{m}-dimensional array with \eqn{1} representing level in the trend component.
#' @param Season  A \eqn{m}-dimensional array  indicating the seasonality for each target series, such as c(12,0).
#' @param vrho  A \eqn{m}-dimensional array of the decay value parameter of the cycle component for each target series, such as c(0,0.99).
#' @param lambda  A \eqn{m}-dimensional array of the frequence parameter of the cycle component for each target series, such as c(0,pi/100).
#' 
#' @examples
#' ###############Setup###########
#' n<-505 #n: sample size
#' m<-2 #m: dimension of target series
#' 
#' cov<-matrix(c(1.1,0.7,0.7,0.9), nrow=2, ncol=2) #covariance matrix of target series 
#' 
#' ###############Regression component###########
#' #coefficients for predictors
#' beta<-t(matrix(c(2,-1.5,0,4,2.5,0,0,2.5,1.5,-1,-2,0,0,-3,3.5,0.5),nrow=2,ncol=8)) 
#' 
#' set.seed(100)
#predictors
#' X1<-rnorm(n,5,5^2)
#' X4<-rnorm(n,-2,5)
#' X5<-rnorm(n,-5,5^2)
#' X8<-rnorm(n,0,100)
#' X2<-rpois(n, 10)
#' X6<-rpois(n, 15)
#' X7<-rpois(n, 20)
#' X3<-rpois(n, 5)
#' X<-cbind(X1,X2,X3,X4,X5,X6,X7,X8) 
#' 
#' ###############Simulated data################
#' set.seed(100)
#' data=sim_data(X=X, beta=beta, cov, k=c(8,8), mu=c(1,1), rho=c(0.6,0.8), 
#'               Dtilde=c(-1,3), Season=c(100,0), vrho=c(0,0.99), lambda=c(0,pi/100))
#' 
#' 
#' 
#'@references
#'\Qiu2018
#'
#'\Ning2021
#'
#'\Jammalamadaka2019
#' @export

setGeneric(
    "sim_data",
    function (X, beta, cov, k, mu, rho, mean_trend=1, sd_trend=0.5, mean_season=20, sd_season=0.5, mean_cycle=20, sd_cycle=0.5, Dtilde, Season, vrho, lambda)
        standardGeneric("sim_data")
)

#' @rdname sim_data
setMethod(
    "sim_data",
    signature=signature(X="array"),
    definition=function (X, beta, cov, k, mu, rho, mean_trend=1, sd_trend=0.5, mean_season=20, sd_season=0.5, mean_cycle=20, sd_cycle=0.5, Dtilde, Season, vrho, lambda)
 {
        tryCatch(
          sim_data.internal(X, beta, cov, k, mu, rho, mean_trend=1, sd_trend=0.5, mean_season=20, sd_season=0.5, mean_cycle=20, sd_cycle=0.5, Dtilde, Season, vrho, lambda)
        )
        
    }
)

sim_data.internal <-
  function(X, beta, cov, k, mu, rho, mean_trend=1, sd_trend=0.5, mean_season=20, sd_season=0.5, mean_cycle=20, sd_cycle=0.5, Dtilde, Season, vrho, lambda){ 
    
    n=dim(X)[1]
    m=dim(beta)[2]
    ###############Trend component###########
    trend=matrix(0,n,m) #Trend component
    delta=matrix(0,n,m) #Slope
    
    for (j in 1:m){
      for(i in 2:n){
        if (mu[j]==T){
          trend[i,j]<-trend[i-1,j]+delta[i-1,j]+rnorm(n=1,mean=mean_trend,sd=sd_trend)
          if(rho[j]!=0){
            delta[i,j]<-Dtilde[j]+rho[j]*(delta[i-1,j]-Dtilde[j])+rnorm(n=1,mean=mean_trend,sd=sd_trend)
          }
        }
      }
    }
    
    
    ###############Seasonal component###########
    sl=matrix(0,n,m) #Seasonal component
    
    for (j in 1:m){
      for(i in 1:n){
        if (Season[j]!=0){
          if(Season[j]>i){
            sl[i,j]<-rnorm(n=1,mean=mean_season,sd=sd_season)
          } else{
            sl[i,j]<- -sum(sl[(i-Season[j]+1):(i-1),j])+rnorm(n=1,mean=mean_season,sd=sd_season)
          }
        }
      }
    }
    
    
    ###############Cycle component###########
    w=matrix(0,n,m); ws=matrix(0,n,m) #Coupled cycle component
    
    for (j in 1:m){
      for(i in 2:n){
        if (vrho[j]!=0){
          w[i,j]<-vrho[j]*cos(lambda[j])*w[i-1,j]+
            vrho[j]*sin(lambda[j])*ws[i-1,j]+rnorm(n=1,mean=mean_cycle,sd=sd_cycle)
          ws[i,j]<- -vrho[j]*sin(lambda[j])*w[i-1,j]+
            vrho[j]*cos(lambda[j])*ws[i-1,j]+rnorm(n=1,mean=mean_cycle,sd=sd_cycle)
        }
      }
    }
    
    
    #regression componenet
    reg<-X%*%beta 
    
    
    ###############Error term###########
    err<-mvrnorm(n=n,mu=rep(0,m),Sigma=cov)  
    
    
    ###############Target series###########
    Y=reg+trend+sl+w+err 
    colnames(Y)<-paste("Y", 1:m, sep = "")
    
    
    ###############Simulated data###########
    simdata = cbind(Y,X)  
    
    return(simdata)
  } 
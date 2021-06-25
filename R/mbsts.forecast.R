#' @name mbsts.forecast
#' @title Specification of time series components
#' @author Jinwen Qiu \email{qjwsnow_ctw@hotmail.com} Ning Ning \email{patricianing@gmail.com}
#' 
#' @description Generate draws from the posterior predictive distribution of a mbsts object. Samples from the posterior predictive distribution of the MBSTS model. 
#' 
#' 
#' @param object An object of the mbsts class created by a call to the mbsts_function function.
#' @param STmodel An object of the SSModel class created by a call to the tsc.setting function.
#' @param newdata A vector or matrix containing the predictor variables to use in making the prediction. This is only required if the mbsts model has a regression component.
#' @param steps An integer value describing the number of time steps ahead to be forecasted. If it is greater than the number of new observations in the newdata, zero values will fill in missing new observations.

#' @return An object of predicted values which is a list containing the following:
#' \item{pred.dist}{An array of draws from the posterior predictive distribution. The first dimension in the array represents time, the second dimension denotes each target series, and the third dimension indicates each MCMC draw.}
#' \item{pred.mean}{A matrix giving the posterior mean of the prediction for each target series.
#' }
#' 

#'@references
#'\Qiu2018
#'
#'\Ning2021
#'
#'\Jammalamadaka2019


#' @export
mbsts.forecast<-
    function(object,STmodel,newdata,steps=1){
        ##extract values from mbsts object
        if(!is.null(STmodel)){
            States<-object@States[dim(object@States)[1],,]
            st.sig2<-object@st.sig2
        }
        if(!is.null(newdata)){
            B.hat<-object@B.hat
        }
        ob.sig2<-object@ob.sig2
        
        if(!is.null(STmodel)){
            ##Initialization of time series components
            ls<-array(0,c(steps,dim(ob.sig2)[1],dim(object@States)[3]))
            st.err<-array(0,c(dim(st.sig2)[1],steps,dim(object@States)[3]))
            ####Predict time series components
            for(i in 1:dim(object@States)[3]){
                for(j in 1:steps){
                    if(j==1){
                        st.err[,j,i]<-t(mvrnorm(n=1,mu=c(rep(0,dim(st.sig2)[1])),
                                                Sigma=diag(st.sig2[,i],dim(st.sig2)[1])))
                        newstate<-STmodel$T[,,1]%*%States[,i]+
                            STmodel$R[,,1]%*%st.err[,j,i]
                    } else{
                        st.err[,j,i]<-t(mvrnorm(n=1,mu=c(rep(0,dim(st.sig2)[1])),
                                                Sigma=diag(st.sig2[,i],dim(st.sig2)[1])))
                        newstate<-STmodel$T[,,1]%*%newstate+
                            STmodel$R[,,1]%*%st.err[,j,i]
                    }
                    ls[j,,i]<-t(STmodel$Z[,,1]%*%newstate)
                }
            }
        }
        
        if(!is.null(newdata)){
            ###check if step is greater than number of observations for newdata
            if(dim(newdata)[1]>=steps){
                newdata<-newdata[1:steps,]
            } else{
                newdata<-rbind(newdata,matrix(0,steps-dim(newdata)[1],dim(newdata)[2]))
            }
            ##Initialization of regression components
            reg<-array(0,c(steps,dim(ob.sig2)[1],dim(B.hat)[3]))
            ##Predict regression components
            for (i in 1:dim(B.hat)[3]){
                reg[,,i]<- newdata%*%B.hat[,,i]
            }
        }
        ####Observation errors
        err<-array(0,c(steps,dim(ob.sig2)[1],dim(ob.sig2)[3]))
        for(j in 1:steps){
            for(i in 1:dim(ob.sig2)[3]){
                err[j,,i]=t(mvrnorm(n=1,mu=c(rep(0,dim(ob.sig2)[1])),
                                    Sigma=ob.sig2[,,i]))  
            }
        }
        if(!is.null(STmodel) & !is.null(newdata)){
            pred.distribution<-ls+reg+err
        } else{
            if(!is.null(STmodel)){
                pred.distribution<-ls+err
            } else{
                pred.distribution<-reg+err
            }
        }
        pred.mean<-apply(pred.distribution,c(1,2),mean)
        
        return(list(pred.dist=pred.distribution,pred.mean=pred.mean))
    }

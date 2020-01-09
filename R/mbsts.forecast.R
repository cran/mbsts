mbsts.forecast <-
function(mbsts,STmodel=NULL,newdata=NULL,steps=1){
    ##extract values from mbsts object
    if(!is.null(STmodel)){
        States<-mbsts$States[dim(mbsts$States)[1],,]
        st.sig2<-mbsts$st.sig2
    }
    if(!is.null(newdata)){
        B.hat<-mbsts$B.hat
    }
    ob.sig2<-mbsts$ob.sig2
    
    if(!is.null(STmodel)){
        ##Initialization of time series components
        ls<-array(0,c(steps,dim(ob.sig2)[1],dim(mbsts$States)[3]))
        st.err<-array(0,c(dim(st.sig2)[1],steps,dim(mbsts$States)[3]))
        ####Predict time series components
        for(i in 1:dim(mbsts$States)[3]){
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
    
    return(list(pred.distribution=pred.distribution,pred.mean=pred.mean))
}

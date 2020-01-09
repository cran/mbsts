mbsts <-
function(Y,X.star=NULL,STmodel=NULL,
                ki=NULL,pii=NULL,
                b=NULL,kapp=0.1,
                R2=0.8,v0=NULL,
                v=0.01,ss=0.01,
                mc=500,burn=50){   
    n=nrow(Y)       #number of observations
    m=ncol(Y)       #number of response variables
    if(!is.null(X.star)){
        I=c(rep(1,dim(X.star)[2]))    #Initialization of indicator function
    }
    
    #Initialization of Sigma_eplison
    if(is.null(STmodel)){
        Sigma2<-diag(0.1,m)
    } else {
        Sigma2<-matrix(STmodel$H,m,m)
    }
       
    ######Initialization of output results
    #####time series components
    if(!is.null(STmodel)){
        States<-array(0,c(dim(Y)[1],dim(STmodel$R)[1],mc-burn))
        st.sig2<-matrix(0,dim(STmodel$R)[2],mc-burn)
    }
    ###regression component
    if(!is.null(X.star)){
        Ind<- matrix(0,dim(X.star)[2],mc-burn)
        beta.hat<- matrix(0,dim(X.star)[2],mc-burn)
        B.hat<-array(0,c(dim(X.star)[2],m,mc-burn))
    }
    ob.sig2<-array(0,c(m,m,mc-burn))  ##Sigma_eplison
    
    for(jj in 1:mc){
        if(!is.null(STmodel)){  #check existence of time series components
            ### draw latent state#### 
            LS<- simulateSSM(STmodel, type = "states")    #draw latent states by smoother
            LS<-matrix(LS,ncol=dim(LS)[2])
            State.error<- simulateSSM(STmodel,type = "eta")  #draw state disturbance
            Obs.error<- matrix(simulateSSM(STmodel,type = "epsilon"),ncol=m)  #draw residuals
            State.error<-matrix(State.error,ncol=dim(State.error)[2])
            
            ##### Draw state component parameters ######
            v.post <-  v+n/2     #posterior degree freedom
            ss.post<-vector()
            State.sigma2<-vector()
            for (i in 1:ncol(State.error)) {
                ss.post[i] <- ss+crossprod(State.error[,i])/2       #posterior sum of square
                State.sigma2[i]<-rigamma(1,alpha=v.post,beta=ss.post[i]) #draw state parameters
            }
        }
        
        if(!is.null(X.star)){ #check existence of regression components
            ####SSVS for drawing gamma, sigma_eplison, beta######
            ##transform design matrix to a larger matrix for computation convenience
            K=ncol(X.star)         #number of predictors
            for (i in 1:K){
                X.star[,i]<-X.star[,i]-mean(X.star[,i])  #demean predictors
            }
            Xtemp1<-X.star[,1:ki[1]]
            Xtemp2<-X.star[,(ki[1]+1):ki[2]]
            X<-as.matrix(bdiag(Xtemp1,Xtemp2))
            if(m>2){
                for (j in 3:m){   
                    Xtemp<-X.star[,(ki[j-1]+1):ki[j]]
                    X<-as.matrix(bdiag(X,Xtemp)) 
                }
            }
            #####Transformation end##############
            if(is.null(STmodel)){
               Y.star<- Y 
            } else {
               Y.star<- Y-LS%*%t(matrix(STmodel$Z,nrow=m))  
               #substract time series component from target series
            }
            for (i in 1:m){
                Y.star[,i]<-Y.star[,i]-mean(Y.star[,i])  #demean response variable
            }
            Y.tilde<-matrix(Y.star,ncol = 1)   #transform to vector form
            
            #transformed system with uncorrelated errors
            U<-chol(Sigma2)    #Cholesky decomposition for observation variance parameter
            UI<-t(solve(U))%x%diag(n)
            X.hat<-UI%*%X            #transformed X
            Y.hat<-UI%*%Y.tilde      #trasnformed y
            
            ###Prior parameters setting
            V0=(v0-m-1)*(1-R2)*var(Y.star)    #prior scale matrix 
            A<-kapp*(t(X)%*%X)/n   #prior information matrix  
            
            ###All posterior parameters 
            V<-t(X.hat)%*%X.hat+A   #posterior co-variance matrix for beta
            N=n+v0       #posterior degree freedom
            gama<-I     #assign value to temporary gamma
            Xindex<-seq(1,K)  #column index for each predictor in X
            
            #####draw each gamma ######
            for(j in sample(Xindex))
            {   
                zero<-0
                ####To make sure at least one predictor selected
                if(sum(gama[-j])==0){
                    zero<-1            #except jth predictor, other predictors not selected
                    Index<-sample(Xindex[-j],1)
                    gama[Index]<-1     #randomly choose one predictor except jth predictor
                }
                
                p<-vector()
                for(value in 1:0)
                { 
                    gama[j]<-value
                    p.gamma<- prod(pii^gama)*prod((1-pii)^(1-gama))     
                    #prior probability for gamma 
                    b.gamma<-b[which(gama==1),]   
                    #prior coefficients for selected predictors
                    X.gamma<-X.hat[,which(gama==1)] #design matrix with selected predictors
                    ####prior parameters
                    A.gamma<-A[which(gama==1),which(gama==1)] 
                    #prior information matrix with selected predictors
                    
                    ####posterior parameters
                    V.gamma<-V[which(gama==1),which(gama==1)]   
                    #posterior co-variance matrix for beta with selected predictors 
                    Z.gamma<-t(X.gamma)%*%Y.hat+A.gamma%*%b.gamma
                    
                    exp.par<-as.vector(t(b.gamma)%*%A.gamma%*%b.gamma-
                                           t(Z.gamma)%*%solve(V.gamma)%*%Z.gamma)
                    ##Set constant value to shrink value for exponent##
                    if(value==1){
                        temp<-exp.par
                    }
                    ###posterior pmf for gamma
                    p[value+1]<-(det(as.matrix(A.gamma))/det(as.matrix(V.gamma)))^(1/2)*
                        p.gamma*exp(-0.5*(exp.par-temp))  
                }
                ##debug due to computation round off
                if (p[2]/(p[1]+p[2])>1){
                    gama[j]<-rbinom(1, 1, prob=1)     
                } else {
                    gama[j]<-rbinom(1, 1, prob=p[2]/(p[1]+p[2]))    #update the gamma[j]
                }
                if(zero==1){
                    gama[Index]<-0
                }
            }
            I<-gama    #assign updated gamma to indicator function
            
            ###All predictors not selected
            if(sum(I)==0){
                Obs.sigma2<-cov(Obs.error)     #assign value to sigma^2_eplison
                beta<-matrix(0,nrow=K,ncol=1)  #assign zero to all coefficients
                B<-matrix(0,nrow=K,ncol=m)     
            }  else{     ###At least one predictor selected
                ######draw sigma and beta############
                b.gamma<-b[which(I==1),]    
                #prior coefficients for selected predictors
                X.gamma<-X.hat[,which(I==1)]   #design matrix with selected predictors
                Xstar.gamma<-X.star[,which(I==1)]   #design matrix with selected predictors
                ####prior parameters
                A.gamma<-A[which(I==1),which(I==1)] 
                #prior information matrix with selected predictors 
                ####posterior parameters
                V.gamma<-V[which(I==1),which(I==1)]   
                #posterior co-variance matrix for beta with selected predictors 
                beta.tilde<-solve(V.gamma)%*%(t(X.gamma)%*%Y.hat+A.gamma%*%b.gamma)
                
                beta<-matrix(0,nrow=K,ncol=1)     #initialization of beta
                beta.temp<-mvrnorm(n=1,mu=beta.tilde,Sigma=solve(V.gamma))    
                #draw beta from multivariate normal distribution
                beta[which(I==1),]<-beta.temp
                
                ####transfrom beta to matrix form
                betai1<-beta[1:ki[1],]
                betai2<-beta[(ki[1]+1):ki[2],]
                B<-as.matrix(bdiag(betai1,betai2))
                if(m>2){
                    for (j in 3:m){   
                        betai<-beta[(ki[j-1]+1):ki[j],]
                        B<-as.matrix(bdiag(B,betai)) 
                    }
                }
                B.gamma<-matrix(B[which(I==1),],ncol=m)  #estimated coefficients with selected predictors
                E.gamma<-Y.star-Xstar.gamma%*%B.gamma   #observed error terms
                SS.post=t(E.gamma)%*%E.gamma+V0    #posterior sum of squares matrix
                Obs.sigma2<-riwish(N,SS.post)  #draw observation variance parameter   
            }
        }
        
        ####assign values to output results during each iteration
        if(jj>burn){
            if(!is.null(STmodel)){
                States[,,(jj-burn)]<-LS
                st.sig2[,(jj-burn)]<-State.sigma2
                if(is.null(X.star)){
                    Obs.sigma2<-cov(Obs.error)     #assign value to sigma^2_eplison
                    ob.sig2[,,jj-burn]<-Obs.sigma2
                }
            }
            if(!is.null(X.star)){
                Ind[,(jj-burn)]<- matrix(I,nrow=dim(X.star)[2])
                beta.hat[,(jj-burn)]<- matrix(beta,nrow=dim(X.star)[2])
                B.hat[,,(jj-burn)]<-B
                ob.sig2[,,jj-burn]<-Obs.sigma2
            }
        }
        
        #####updating state space model parameters
        if(!is.null(STmodel)){
            if(!is.null(X.star)){
                rowname<-row.names(STmodel$y)
                STmodel$y[]<- ts(Y-X.star%*%B,names = colnames(STmodel$y)) 
                #target series
                row.names(STmodel$y)<-rowname
                ####variance-covariance matrix for observation errors
                Sigma2 <- matrix(Obs.sigma2,m,m)
                STmodel$H[]<-array(Sigma2,c(m,m,1))
            } else {
                ####variance-covariance matrix for observation errors
                Sigma2 <- matrix(cov(Obs.error),m,m)
                STmodel$H[]<-array(Sigma2,c(m,m,1))
            }
            ##variance-covariance matrix for time series components
            STmodel$Q[]<- array(diag(State.sigma2,nrow = dim(STmodel$R)[2])
                                ,c(dim(STmodel$R)[2],dim(STmodel$R)[2],1))
        } else{
            Sigma2 <- matrix(Obs.sigma2,m,m)
        }
    }
    
    ##return output results
    if(!is.null(X.star) & !is.null(STmodel)){
        return (list(Ind=Ind, beta.hat=beta.hat, B.hat=B.hat,
            ob.sig2=ob.sig2, States=States, st.sig2=st.sig2))
    } else {
        if(!is.null(STmodel)){
            return (list(ob.sig2=ob.sig2, States=States, st.sig2=st.sig2))
        } else {
            return (list(Ind=Ind,beta.hat=beta.hat,B.hat=B.hat,ob.sig2=ob.sig2))
        }
    }
}

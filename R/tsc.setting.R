tsc.setting <-
function(Y,mu=NULL,rho=NULL,S=NULL,vrho=NULL,lambda=NULL){
    ####Space state model setting
    ### y_t=Z*alpha_t+eplsion_t   alpha_t+1=T_t*alpha_t+R_t*alta_t
    ### eplison_t~N(0,H_t)   alta_t~N(0,Q_t)
    ### alpha_1~N(a_1,P_1) 
    
    ###Z:output matrix
    ###T_t:transition matrix
    ###R_t:control matrix
    
    m=dim(Y)[2] # number of target series
    
    # number of latent states with error
    ms1=length(which(mu==T))+length(which(rho!=0))+length(which(S!=0))+
        2*length(which(vrho!=0))
    # number of latent states
    ms2=length(which(mu==T))+length(which(rho!=0))+
        (sum(S)-length(which(S!=0)))+2*length(which(vrho!=0))
    
    ######initilization of matrices
    temp<-vector()
    output<-vector()
    transition<-vector()
    control<-vector()
    
    #######assign values to these matrices
    ####trend#############
    for (i in 1:m){
         if(!is.null(rho)){  ###trend slope
             if(rho[i]!=0){
                 temp<-c(temp,c(1,0))
                 
                 if(length(transition)==0){
                     transition=matrix(c(1,0,1,rho[i]),2,2)
                 } else{
                     transition<-as.matrix(bdiag(transition,
                                                 matrix(c(1,0,1,rho[i]),2,2)))
                 }
                 
                 if(length(control)==0){
                     control=diag(1,2)
                 } else{
                     control<-as.matrix(bdiag(control,diag(1,2)))
                 }
             } else {
                 if(mu[i]==T){
                     temp<-c(temp,1)
                     
                     if(length(transition)==0){
                         transition=matrix(1)
                     } else{
                         transition<-as.matrix(bdiag(transition,1))
                     }
                     
                     if(length(control)==0){
                         control=matrix(1)
                     } else{
                         control<-as.matrix(bdiag(control,1))
                     }
                 }
             }
         } else{  
             if(!is.null(mu)){  #local level 
                 if(mu==T){
                     temp<-c(temp,1)
                     
                     if(length(transition)==0){
                         transition=matrix(1)
                     } else{
                         transition<-as.matrix(bdiag(transition,1))
                     }
                     
                     if(length(control)==0){
                         control=matrix(1)
                     } else{
                         control<-as.matrix(bdiag(control,1))
                     }
                 }
             }
           }
         #########seasonal######
         if (!is.null(S)){
             if (S[i]!=0){
                 temp<-c(temp,c(1,rep(0,S[i]-2)))
                 
                 if(length(transition)==0){
                     transition=rbind(matrix(c(rep(-1,S[i]-1)),1,S[i]-1),
                            matrix(c(diag(1,S[i]-2),rep(0,S[i]-2)),S[i]-2,S[i]-1))
                 } else{
                     transition<-as.matrix(bdiag(transition,
                             rbind(matrix(c(rep(-1,S[i]-1)),1,S[i]-1),
                             matrix(c(diag(1,S[i]-2),rep(0,S[i]-2)),S[i]-2,S[i]-1))))
                 }
                 
                 if(length(control)==0){
                     control=as.matrix(c(1,rep(0,S[i]-2)))
                 } else{
                     control<-as.matrix(bdiag(control,c(1,rep(0,S[i]-2))))
                 }
            }
         }
         ######cycle####
         if(!is.null(vrho)){
             if (vrho[i]!=0){
                 temp<-c(temp,c(1,0))
                 
                 if(length(transition)==0){
                    transition=matrix(c(vrho[i]*cos(lambda[i]),vrho[i]*sin(lambda[i]),
                                -vrho[i]*sin(lambda[i]),vrho[i]*cos(lambda[i])),
                                2,2,byrow = TRUE) 
                 } else{
                     transition<-as.matrix(bdiag(transition,
                            matrix(c(vrho[i]*cos(lambda[i]),vrho[i]*sin(lambda[i]),
                            -vrho[i]*sin(lambda[i]),vrho[i]*cos(lambda[i])),
                            2,2,byrow = TRUE)))
                 }
                 
                 if(length(control)==0){
                     control=matrix(c(1,0,0,1),2,2)
                 } else{
                     control<-as.matrix(bdiag(control,matrix(c(1,0,0,1),2,2)))
                 }
             }
         }
         
         ####construct output matrix######
         if(length(output)==0){
             output<-t(temp)
             temp<-vector()
         } else{
             if(length(temp)!=0){
                 output<-as.matrix(bdiag(output,t(temp)))
                 temp<-vector() 
             }
         }
    }
    #######variance-covariance matrices for errors
    Qt <-diag(0.01,ms1)
    a1 <- matrix(rep(0,ms2),ms2,1)
    P1 <- matrix(0,ms2,ms2)
    P1inf <- diag(ms2)
    Ht=Ht <- diag(0.1,m)   #variance-covariance matrices for observation errors
    
    #####Customized Structural time series model###
    STmodel <- SSModel(Y~ -1+SSMcustom(Z = output, T = transition, R = control, 
               Q = Qt, a1 = a1, P1 = P1,P1inf = P1inf), H = Ht)
    ########
    return(STmodel)
}

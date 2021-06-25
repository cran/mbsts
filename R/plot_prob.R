
#' @title Plot Inclusion Probabilities
#' @author Jinwen Qiu \email{qjwsnow_ctw@hotmail.com} Ning Ning \email{patricianing@gmail.com}
#' @description Plots of the empirical inclusion probabilities for predictors of each target series, based on a user-defined threshold probability. For example, one predictor is selected \eqn{100} times in \eqn{200} MCMC draws (after discard burn-in draws), the empirical inclusion probability for that predictor is \eqn{0.5}. If the user-defined threshold probability less than or equal to \eqn{0.5}, then this predictor will show in the plot.
#' 
#' @importFrom ggplot2 ggplot aes geom_bar scale_fill_manual geom_text position_stack labs theme element_text
#' @docType methods
#' 
#' @param object An object of the mbsts class created by a call to the mbsts_function function.
#' @param title NULL or A character vector whose entries are titles for the inclusion probability plots generated for each target series, such as c("Inclusion Probabilities for y1", "Inclusion Probabilities for y2"). If Null, the output is c("y1","y2",...).
#' @param prob.threshold A numerical value used as the threshold to only include predictors whose inclusion probabilities are higher than it in the plot. The default value is \eqn{0.2}.
#' @param varnames NULL or A character vector whose entries are the variable names for predictors, such as c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8"). If Null, the output is c("x11","x12",...,"x21","x22",...). 

#'@references
#'\Qiu2018
#'
#'\Ning2021
#'
#'\Jammalamadaka2019

#' @export
setGeneric(
    "plot_prob",
    function (object,title=NULL,prob.threshold=0.2,
              varnames=NULL)
        standardGeneric("plot_prob")
)

#' @rdname plot_prob
setMethod(
    "plot_prob",
    signature=signature(object="mbsts"),
    definition=function (
        object,
        title=NULL,prob.threshold=0.2,
        varnames=NULL) {
        tryCatch(
            plot_prob.internal(
                object,
                title=title,prob.threshold=prob.threshold,
                varnames=varnames
            )
        )
        
    }
)

plot_prob.internal <-
    function(object,
             title=NULL,prob.threshold=0.2,
             varnames=NULL){
        
        Indicator=object@Ind
        coef=object@beta.hat
        seq=object@ki
            
        allplots <- vector('list',length(seq))
        if (is.null(title)){
            for (i in 1:length(seq)){
                title=c(title,paste0('y',toString(i)))
            }
        }
        
        if (is.null(varnames)){
            for (i in 1:length(seq)) {
                if (i==1){
                    for (j in 1:seq[i]){
                        varnames=c(varnames,paste0('x',toString(i),toString(j)))
                    }
                }  else{
                    for (j in 1:(seq[i]-seq[i-1])){
                        varnames=c(varnames,paste0('x',toString(i),toString(j)))
                    }
                }
            } 
        }
        
        
        for (i in 1:length(seq)){
            
            if (i==1){
                
                name<-data.frame(varnames[1:seq[i]])
                prob<-round(apply(Indicator[1:seq[i],],1,mean),2)
                coef.ave <- apply(coef[1:seq[i],],1,mean)
                
            } else{
                
                name<-data.frame(varnames[(seq[i-1]+1):seq[i]])
                prob<-round(apply(Indicator[(seq[i-1]+1):seq[i],],1,mean),2)
                coef.ave <- apply(coef[(seq[i-1]+1):seq[i],],1,mean)
            }
            
            sign<-data.frame(factor(ifelse(coef.ave>0,"positive","negative")))
            plot.data= cbind(name,prob,coef.ave,sign)
            colnames(plot.data)<-c("name","prob","coef","sign")
            
            allplots[[i]]<-local({
                p <- ggplot(plot.data[which(plot.data$prob>=prob.threshold),],
                            aes( x= name,y= prob,fill=sign ,label = prob))+
                    geom_bar(stat = "identity")+coord_flip()+
                    scale_fill_manual(values  = c("blue", "red"))+
                    geom_text(size = 3.5, position = position_stack(vjust = 1.08))+
                    labs(title=title[i], 
                         y="Probability",
                         x="Variable Name",color=NULL)+
                    theme(text = element_text(size=12),
                          axis.text = element_text(size=14),
                          plot.title = element_text(hjust = 0.4,size=16),
                          legend.text=element_text(size=12),
                          legend.title=element_text(size=12))
            })
            
        }
       
        return (allplots)
    } 


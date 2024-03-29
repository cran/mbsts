
#' @title Plot Posterior State Components
#' @author Jinwen Qiu \email{qjwsnow_ctw@hotmail.com} Ning Ning \email{patricianing@gmail.com}
#' @description Plots of the mean of posterior state components of each target series, which is generated by the model training procedure of the MBSTS model.
#' 
#' @importFrom ggplot2 ggplot aes geom_line scale_colour_discrete facet_wrap labs theme
#' @importFrom BBmisc normalize
#' @importFrom reshape2 melt
#' @docType methods
#' 
#' @param object An object of the mbsts class created by a call to the mbsts_function function.
#' @param slope A logical vector indicating whether there is trend for each target series, such as c(T,T).
#' @param local A logical vector indicating whether there is local level for each target series, such as c(T,T).
#' @param season A numerical vector indicating the seasonality for each target series, such as c(12,0).
#' @param cyc A logical vector indicating whether there is a cycle component for each target series, such as c(F,T).
#' @param time Null or a data frame for time index of the time series. The default value is data.frame(seq(1,n)).
#' @param title NULL or a character vector whose entries are titles for the plots of target series' posterior state components, such as c("Posterior State Components of y1", "Posterior State Components of y2"). The default is c("y1","y2",...).
#' @param component_selection A character variable whose value must be one of "All", "Trend", "Seasonal", "Cycle", and "Regression". Here, "Trend" means the trend component only and "All" means all the components.

#'@references
#'\Qiu2018
#'
#'\Ning2021
#'
#'\Jammalamadaka2019


#' @export
setGeneric(
  "plot_comp",
  function (object,slope,local,season,cyc,time=NULL, 
            title=NULL,component_selection="All")
    standardGeneric("plot_comp")
)

#' @rdname plot_comp
setMethod(
  "plot_comp",
  signature=signature(object="mbsts"),
  definition=function (
    object,
    slope,local,season,cyc,time=NULL, 
    title=NULL,component_selection="All"
    ) {
    tryCatch(
      plot_comp.internal(
        object,
        slope=slope,local=local,season=season,cyc=cyc,
        time=time,title=title,component_selection=component_selection
      )
    )
  }
)


################ plot state component ##############
plot_comp.internal <-
  function(object,
           slope,local,season,cyc,time=NULL, 
           title=NULL,component_selection="All"){ 
    
    X=object@Xtrain
    states=object@States
    Bhat=object@B.hat
    n=object@ntrain
    m=object@mtrain
    
    allplots <- vector('list',m)
    indplots <- vector('list',m)
    
    if (is.null(title)){
      for (i in 1:m){
        title=c(title,paste0('y',toString(i)))
      }
    }
    
    if (!is.null(local)|!is.null(season)|!is.null(cyc)){
      state.ave<-apply(states,c(1,2),mean)
    }
    
    
    if (!is.null(X)){
      B.ave <- apply(Bhat,c(1,2),mean)
      reg.all<-data.frame(X%*%B.ave)
    }
    
    pos=1
    
    for (i in 1:m){
      #i=2
      plot.data<-data.frame(seq(1,n))
      
      if(slope[i]==T){
        trend<- data.frame(state.ave[,pos])
        trend<-normalize(trend,method="range",range=c(-1,1))
        plot.data<-cbind(plot.data,trend)
        pos=pos+2
      } else {
        if (local[i]==T){
          trend<-data.frame(state.ave[,pos])
          trend<-normalize(trend,method="range",range=c(-1,1))
          plot.data<-cbind(plot.data,trend)
          pos=pos+1
        } else {
          trend=NULL
        }
      }
      
      if(season[i]!=0){
        seasonal<-data.frame(state.ave[,pos])
        seasonal<-normalize(seasonal,method="range",range=c(-1,1))
        plot.data<-cbind(plot.data,seasonal)
        pos = pos+season[i]-1
      } else {
        seasonal=NULL
      } 
      
      if(cyc[i]==T){
        cycle<-data.frame(state.ave[,pos])
        cycle<-normalize(cycle,method="range",range=c(-1,1))
        plot.data<-cbind(plot.data,cycle)
        pos = pos+2
      } else {
        cycle=NULL
      }                
      
      if (!is.null(X)){
        reg = data.frame(reg.all[,i])
        reg<-normalize(reg,method="range",range=c(-1,1))
        plot.data<-cbind(plot.data,reg)
      } else {
        reg=NULL
      }
      
      
      cname=c("time")
      if(!is.null(trend)){
        cname<-c(cname,"Trend")
      }
      
      if(!is.null(seasonal)){
        cname<-c(cname,"Seasonal")
      }
      
      if(!is.null(cycle)){
        cname<-c(cname,"Cycle")
      }
      
      if(!is.null(reg)){
        cname<-c(cname,"Regression")
      }
      
      colnames(plot.data)<-cname
      row.names(plot.data)<-as.character(seq(1,n))
      plot.data2<-melt(plot.data,id=c("time"))
      
      if (component_selection=="All") {
        value<-plot.data2[,3]
        variable<-plot.data2[,2]
        allplots[[i]]<-ggplot(plot.data2, aes(x = time,y=value,colour=variable)) +
          geom_line(aes(colour=variable)) +
          scale_colour_discrete(guide = "none") +
          facet_wrap(~variable,ncol= 1)+
          labs(title=title[i], x="Time",
               y="State Components", color=NULL)+
          theme(text = element_text(size=13),
                axis.text = element_text(size=14),
                plot.title = element_text(hjust = 0.4,size=14))
      }else {
        a=which(plot.data2[2]==component_selection)
        plot.data2.ind=plot.data2[a,]
        if(dim(plot.data2.ind)[1]==0){
          indplots[[i]]<-NULL
        }else{
          value<-plot.data2.ind[,3]
          variable<-plot.data2.ind[,2]
          indplots[[i]]<-ggplot(plot.data2.ind, aes(x = time,y=value,colour=variable)) +
            geom_line(aes(colour=variable)) +
            scale_colour_discrete(guide = "none") +
            facet_wrap(~variable,ncol= 1)+
            labs(title=title[i], x="Time",
                 y="State Component", color=NULL)+
            theme(text = element_text(size=13),
                  axis.text = element_text(size=14),
                  plot.title = element_text(hjust = 0.4,size=15))
          
        }
      }
      
    }
    
    if (component_selection=="All") {
      return (allplots)
    } else {
      return (indplots)
    }
    
  } 





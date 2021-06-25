
#' @title Regression parameter estimation by the MBSTS Model
#' @author Jinwen Qiu \email{qjwsnow_ctw@hotmail.com} Ning Ning \email{patricianing@gmail.com}
#' @description Generate feature selection and parameter estimation results of a mbsts object. Provide means and standard deviations of parameter estimation results for selected features.
#' 
#' @docType methods
#' @importFrom matrixStats rowSds
#' 
#' @param object An object of the mbsts class created by a call to the mbsts_function function.
#' @param prob.threshold A numerical value used as the threshold to only include predictors whose inclusion probabilities are higher than it in the plot. The default is \eqn{0.2}.#' @param prob.threshold A numerical value used as the threshold to only include predictors whose inclusion probabilities are higher than it in the plot. The default value is \eqn{0.2}.


#' @return A list with the following components
#' \item{index}{An array of feature selection results.}
#' \item{para.est.mean}{An array of means of parameter estimation values of selected features.}
#' \item{para.est.sd}{An array of standard deviations of parameter estimation values of selected features.}


#'@references
#'\Qiu2018
#'
#'\Ning2021
#'
#'\Jammalamadaka2019


#' @export
setGeneric(
  "para.est",
  function (object,prob.threshold=0.2)
    standardGeneric("para.est")
)

#' @rdname para.est
setMethod(
  "para.est",
  signature=signature(object="mbsts"),
  definition=function (
    object,
    prob.threshold=0.2
  ) {
    tryCatch(
      para.est.internal(
        object,
        prob.threshold=prob.threshold
      )
    )
  }
)

para.est.internal <-function(object,prob.threshold=0.2){

index=which(rowMeans(object@Ind)>prob.threshold)
para.est.mean=round(rowMeans(object@beta.hat)[index],digits = 4)
para.est.sd=round(rowSds(object@beta.hat)[index],digits = 4)

return (list(index=index, para.est.mean=para.est.mean, para.est.sd=para.est.sd))
}
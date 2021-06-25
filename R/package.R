## package description

#' @docType package
#' @name mbsts-package
#' @title Multivariate Bayesian Structural Time Series
#'
#' @description Tools for data analysis with multivariate Bayesian structural time series (MBSTS) models.  Specifically, the package provides facilities for implementing general structural time series models, flexibly adding on different time series components (trend, season, cycle, and regression), simulating them, fitting them to multivariate correlated time series data, conducting feature selection on the regression component.
#'
#' @section Documentation:
#' \pkg{mbsts} is described in Ning and Qiu (2021).
#'
#' @section License:
#' \pkg{mbsts} is provided under the \acronym{LGPL-2.1} License.
#'
#' @references  
#' \Qiu2018
#' 
#' \Jammalamadaka2019
#' 
#' \Ning2021
#' 
#' @author Jinwen Qiu, Ning Ning
#' @import methods
#' @importFrom Matrix bdiag
#' @importFrom KFAS simulateSSM
#' @importFrom pscl rigamma
#' @importFrom MASS mvrnorm
#' @importFrom MCMCpack riwish
#' @importFrom stats cov rbinom ts var
#' @importFrom matrixStats rowSds
#' @import ggplot2 
#' @importFrom BBmisc normalize
#' @importFrom reshape2 melt
#' @importFrom KFAS SSModel SSMcustom
#' @importFrom Matrix bdiag
NULL        

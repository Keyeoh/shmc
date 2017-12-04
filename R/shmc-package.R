#' Simple hierarchical model of 5hmC
#'
#' @docType package
#' @name shmc
#' @useDynLib shmc, .registration = TRUE
#'
#' @import methods
#' @import Rcpp
#'
NULL

# quiets concerns of R CMD check
if(getRversion() >= '2.15.1')  {
  utils::globalVariables(c(
    'Contrast',
    'ESS',
    'Group',
    'Param',
    'Sd',
    'State',
    'Value'
  ))
}

#' @include covarianceStructure3D.R
NULL

#### covarianceModel3D class ####
#' Covariance model
#'
#' Representation of a spatial continuity model.
#'
#' @slot nugget The nugget effect or noise variance.
#' @slot nugget.dir The nugget effect for directional data.
#' @slot structures A \code{list} with one or more
#' \code{covarianceStructure3D} objects
#' @slot total.var The total variance of the covariance structures,
#' excluding the nugget effect.
#'
#' @seealso \code{\link{covarianceStructure3D-init}},
#' \code{\link{CovarianceMatrix}}
#'
#' @export covarianceModel3D
covarianceModel3D <- setClass(
  "covarianceModel3D",
  slots = c(nugget = "numeric",
            nugget.dir = "numeric",
            structures = "list",
            total.var = "numeric"),
  validity = function(object){
    if (length(nugget) != 1)
      stop("nugget must be a scalar")
    if (length(nugget.dir) != 1)
      stop("nugget.dir must be a scalar")
    if (!all(rapply(structures, class) == "covarianceStructure3D")){
      stop("all structures must be of class 'covarianceStructure3D'")
    }
  }
)

#### initialize ####
#' Covariance model
#'
#' Representation of a spatial continuity model.
#'
#' @param nugget The nugget effect or noise variance.
#' @param nugget.dir The nugget effect for directional data.
#' @param structures Objects of the \code{covarianceStructure3D} class
#' describing the desired covariance structure. A single object or multiple,
#' concatenated objects.
#'
#' @name covarianceModel3D
#'
#' @seealso \code{\link{covarianceStructure3D-class}}
covarianceModel3D <- function(nugget, structures, nugget.dir = 0){
  if (class(structures) != "list") structures <- list(structures)
  total.var <- sum(sapply(.Object@structures, function(s) s@contribution))
  new("covarianceModel3D", nugget = nugget, nugget.dir = nugget.dir,
      structures = structures, total.var = total.var)
}

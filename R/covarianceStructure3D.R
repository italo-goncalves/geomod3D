#' @include generics.R
#' @include utils.R
NULL

#### covarianceStructure3D class ####
#' Covariance structure
#'
#' Representation of one structure in a spatial continuity model.
#'
#' @slot type The type of covariance function represented.
#' @slot contribution The variance associated with this structure.
#' @slot maxrange,midrange,minrange The structure's range in three orthogonal
#' directions, or the three semi-axes of the covariance ellipsoid.
#' @slot azimuth,dip,rake Orientation of the covariance ellipsoid in geological
#' coordinates.
#' @slot power The power parameter for the covariance functions that feature
#' one.
#' @slot anis An anisotropy matrix.
#' @slot covfun,covd1,covd2,varfun,vard1,vard2 Functions to calculate the
#' covariance and variogram matrices between a variable of interest and/or its
#' derivatives.
#'
#' @seealso \code{\link{covarianceStructure3D-init}},
#' \code{\link{CovarianceMatrix}}
#'
#' @export covarianceStructure3D
covarianceStructure3D <- setClass(
  "covarianceStructure3D",
  slots = c(type = "character",
            contribution = "numeric",
            maxrange = "numeric",
            midrange = "numeric",
            minrange = "numeric",
            azimuth = "numeric",
            dip = "numeric",
            rake = "numeric",
            power = "numeric",
            anis = "matrix",      # anisotropy matrix
            covfun = "function",  # value/value covariance
            covd1 = "function",   # value/derivative covariance
            covd2 = "function"),  # derivative/derivative covariance
  validity = function(object) {
    if (!(object@type %in% c("gaussian",
                             "exponential",
                             "spherical",
                             "power",
                             "cubic",
                             "matern1",
                             "matern2",
                             "cauchy")))
      stop(paste0("Structure type '",object@type, "' not supported."))
    if (object@contribution <= 0)
      stop("Contribution must be greater than zero")
    if (min(c(object@maxrange, object@midrange, object@minrange)) <= 0)
      stop("All ranges must be greater than zero")
    if (object@type == "power")
      if (object@power <= 0 || object@power >= 2)
        stop("Power must be in the interval ]0,2[ for power function")
    if (object@type == "cauchy" && object@power <= 0)
      stop("Power must be positive for cauchy function")
    return(TRUE)
  }
)

#### initialize ####
#' Covariance structure
#'
#' Representation of one structure in a spatial continuity model.
#'
#' @param type The type of covariance function represented (see below).
#' @param contribution The variance associated with this structure.
#' @param maxrange,midrange,minrange The structure's range in three orthogonal
#' directions, or the three semi-axes of the covariance ellipsoid.
#' @param azimuth,dip,rake Orientation of the covariance ellipsoid in geological
#' coordinates.
#' @param power The power parameter for the covariance functions that feature
#' one.
#'
#' @section Covariance functions:
#'
#' The following functions are available for the \code{type} parameter. In the
#' equations below, \eqn{d} is the anisotropy-corrected distance between two points
#' and \eqn{p} is the \code{power} parameter. Functions marked with an asterisk
#' support modeling with structural data.
#'
#' \describe{
#'   \item{\code{gaussian}*}{\deqn{\exp{-3d^2}}{exp(-3*d^2)}}
#'   \item{\code{exponential}}{\deqn{\exp{-3d}}{exp(-3*d)}}
#'   \item{\code{spherical}}{\deqn{1 - 1.5 d + 0.5 d^3, 0 \leq d \leq 1}{1 - 1.5*d + 0.5*d^3, 0 <= d <= 1}}
#'   \item{\code{cubic}*}{\deqn{1 - 7d^2 + \frac{35}{4}d^3 - \frac{7}{2}d^5 + \frac{3}{4}d^7, 0 \leq d \leq 1}{1 - 7*d^2 + 35/4*d^3 - 7/2*d^5 + 3/4*d^7, 0 <= d <= 1}}
#'   \item{\code{matern1}*}{\deqn{(1 + 5d)\exp{-5d}}{(1 + 5*d) * exp(-5*d)}}
#'   \item{\code{matern2}*}{\deqn{(1 + 6d + 12d^2)\exp{-6d}}{(1 + 6*d + 12*d^2) * exp(-6*d)}}
#'   \item{\code{cauchy}*}{\deqn{(1 + d^2)^{-p}, p > 0}{(1 + d^2)^(-p), p > 0}}
#' }
#'
#' @name covarianceStructure3D
#'
#' @seealso \code{\link{covarianceStructure3D-class}}
covarianceStructure3D <- function(type = c("gaussian", "exponential",
                                           "spherical", "cubic", "matern1",
                                           "matern2", "cauchy"),
                                  contribution,
                                  maxrange, midrange = maxrange,
                                  minrange = midrange,
                                  azimuth = 0, dip = 0, rake = 0,
                                  power = 1){
  type <- type[1]
  anis <- solve(anisotropy3D(maxrange, midrange, minrange, azimuth, dip, rake))
  if (type == "gaussian"){
    covfun <- function(u, v = u){
      A <- anis
      v1 <- t(A %*% t(v))
      u1 <- t(A %*% t(u))
      d <- .vectorized_pdist(u1, v1)
      contribution * exp(-3*(d^2))
    }
    covd1 <- function(u, v, dir1){
      contribution * .covd1_gaussian(u, v, dir1, anis)
    }
    covd2 <- function(u, v = u, dir1, dir2){
      contribution * .covd2_gaussian(u, v, dir1, dir2, anis)
    }
  }
  else if (type == "exponential"){
    covfun <- function(u, v = u){
      A <- anis
      v1 <- t(A %*% t(v))
      u1 <- t(A %*% t(u))
      d <- .vectorized_pdist(u1, v1)
      contribution * exp(-3 * d)
    }
    covd1 <- function(u, v, dir1){
      stop("exponential function does not support directional data")
    }
    covd2 <- function(u, v = u, dir1, dir2){
      stop("exponential function does not support directional data")
    }
  }
  else if (type == "spherical"){
    covfun <- function(u, v = u){
      A <- anis
      v1 <- t(A %*% t(v))
      u1 <- t(A %*% t(u))
      d <- .vectorized_pdist(u1, v1)
      d[d > 1] <- 1
      contribution * (1 - 1.5 * d + 0.5 * d^3)
    }
    covd1 <- function(u, v, dir1){
      stop("spherical function does not support directional data")
    }
    covd2 <- function(u, v = u, dir1, dir2){
      stop("spherical function does not support directional data")
    }
  }
  else if (type == "cubic"){
    covfun <- function(u, v = u){
      A <- anis
      v1 <- t(A %*% t(v))
      u1 <- t(A %*% t(u))
      d <- .vectorized_pdist(u1, v1)
      d[d > 1] <- 1
      contribution * (1 - 7 * d^2 + 35/4 * d^3 - 7/2 * d^5 + 3/4 * d^7)
    }
    covd1 <- function(u, v, dir1){
      contribution * .covd1_cubic(u, v, dir1, anis)
    }
    covd2 <- function(u, v = u, dir1, dir2){
      contribution * .covd2_cubic(u, v, dir1, dir2, anis)
    }
  }
  else if (type == "matern1"){
    covfun <- function(u, v = u){
      A <- anis
      v1 <- t(A %*% t(v))
      u1 <- t(A %*% t(u))
      d <- .vectorized_pdist(u1, v1)
      contribution * (1 + 5 * d) * exp(-5 * d)
    }
    covd1 <- function(u, v, dir1){
      contribution * .covd1_matern1(u, v, dir1, anis)
    }
    covd2 <- function(u, v = u, dir1, dir2){
      contribution * .covd2_matern1(u, v, dir1, dir2, anis)
    }
  }
  else if (type == "matern2"){
    covfun <- function(u, v = u){
      A <- anis
      v1 <- t(A %*% t(v))
      u1 <- t(A %*% t(u))
      d <- .vectorized_pdist(u1, v1)
      contribution * (1 + 5*d + 12*d^2) * exp(-6 * d)
    }
    covd1 <- function(u, v, dir1){
      contribution * .covd1_matern2(u, v, dir1, anis)
    }
    covd2 <- function(u, v = u, dir1, dir2){
      contribution * .covd2_matern2(u, v, dir1, dir2, anis)
    }
  }
  else if (type == "cauchy"){
    covfun <- function(u, v){
      A <- anis
      v1 <- t(A %*% t(v))
      u1 <- t(A %*% t(u))
      beta <- 0.05^(-1/.Object@power) - 1
      d <- .vectorized_pdist(u1, v1)
      contribution * (1+beta*d^2) ^ (- power)
    }
    covd1 <- function(u, v, dir1){
      contribution * .covd1_cauchy(u, v, dir1, anis, power)
    }
    .Object@covd2 <- function(u, v, dir1, dir2){
      contribution * .covd2_cauchy(u, v, dir1, dir2, anis, power)
    }
  }

  # end
  new("covarianceStructure3D", type = type, contribution = contribution,
      maxrange = maxrange, midrange = midrange, minrange = minrange,
      azimuth = azimuth, dip = dip, rake = rake, power = power,
      anis = anis, covfun = covfun, covd1 = covd1, covd2 = covd2)
}

#### show ####
setMethod(
  f = "show",
  signature = "covarianceStructure3D",
  definition = function(object){
    cat("Object of class ", class(object), "\n\n", sep = "")
    cat("Type:", object@type,"\n")
    cat("Contribution =", object@contribution,"\n")
    cat("Maximum range =", object@maxrange,"\n")
    cat("Intermediate range =", object@midrange,"\n")
    cat("Minimum range =", object@minrange,"\n")
    cat("Orientation: azimuth =", object@azimuth,
        "dip =", object@dip, "rake =", object@rake, "\n")
    if(object@type == "power"){
      cat("Power =", object@power,"\n")
    }
  }
)

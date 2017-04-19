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
#' @details A \code{list} with multiple \code{covarianceStructure3D} objects
#' can pe passed on to functions that work with a covariance model. In this
#' case the contributions of each object are summed up.
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
            covd2 = "function",   # derivative/derivative covariance
            varfun = "function",  # value/value variogram
            vard1 = "function",   # value/derivative variogram
            vard2 = "function"),  # derivative/derivative variogram
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
#' @name covarianceStructure3D-init
#'
#' @seealso \code{\link{covarianceStructure3D-class}}
setMethod(
  f = "initialize",
  signature = "covarianceStructure3D",
  definition = function(.Object,
                        type = c("gaussian", "exponential", "spherical",
                                 "cubic", "matern1", "matern2", "cauchy"),
                        contribution,
                        maxrange, midrange = maxrange,
                        minrange = midrange,
                        azimuth = 0, dip = 0, rake = 0,
                        power = 1){
    type <- type[1]
    .Object@type <- type
    .Object@contribution <- contribution
    .Object@maxrange <- maxrange
    .Object@midrange <- midrange
    .Object@minrange <- minrange
    .Object@azimuth <- azimuth
    .Object@dip <- dip
    .Object@rake <- rake
    .Object@power <- power
    validObject(.Object)
    .Object@anis <- solve(anisotropy3D(.Object@maxrange,
                                       .Object@midrange,
                                       .Object@minrange,
                                       .Object@azimuth,
                                       .Object@dip,
                                       .Object@rake
    ))
    if (.Object@type == "gaussian"){
      .Object@covfun <- function(u, v = u){
        A <- .Object@anis
        v1 <- t(A %*% t(v))
        u1 <- t(A %*% t(u))
        d <- .vectorized_pdist(u1, v1)
        .Object@contribution * exp(-3*(d^2))
      }
      .Object@covd1 <- function(u, v, dir1){
        .Object@contribution *
          .covd1_gaussian(u, v, dir1, .Object@anis)
      }
      .Object@covd2 <- function(u, v = u, dir1, dir2){
        .Object@contribution *
          .covd2_gaussian(u, v, dir1, dir2, .Object@anis)
      }
      .Object@varfun <- function(u, v = u){
        .Object@contribution - .Object@covfun(u, v)
      }
      .Object@vard1 <- function(u, v, dir1){
        - .Object@covd1(u, v, dir1)
      }
      .Object@vard2 <- function(u, v, dir1, dir2){
        - .Object@covd2(u, v, dir1, dir2)
      }
    }
    else if (.Object@type == "exponential"){
      .Object@covfun <- function(u, v = u){
        A <- .Object@anis
        v1 <- t(A %*% t(v))
        u1 <- t(A %*% t(u))
        d <- .vectorized_pdist(u1, v1)
        .Object@contribution * exp(-3*d)
      }
      .Object@covd1 <- function(u, v, dir1){
        eps <- .Machine$double.eps
        A <- .Object@anis
        AtA <- t(A) %*% A
        cont <- .Object@contribution
        u <- apply(u, 1, list)
        u <- lapply(u, unlist)
        v <- apply(v, 1, list)
        v <- lapply(v, unlist)
        dir1 <- apply(dir1, 1, list)
        dir1 <- lapply(dir1, unlist)
        ans <- lapply(u, function(u){
          mapply(function(u, v, dir1){
            h <- A %*% matrix(u - v, 3, 1)
            d <- sqrt(sum(h^2))
            K <- exp(-3*d)
            h <- matrix(u - v, 3, 1)
            grad <- cont * 3 * K /(d+eps) *
              AtA %*% h
            sum(grad * dir1)
          }, list(u), v, dir1)
        })
        matrix(unlist(ans), length(u), length(v),
               byrow = T)
      }
      .Object@covd2 <- function(u, v = u, dir1, dir2){
        eps <- .Machine$double.eps
        A <- .Object@anis
        AtA <- t(A) %*% A
        cont <- .Object@contribution
        u <- apply(u, 1, list)
        u <- lapply(u, unlist)
        v <- apply(v, 1, list)
        v <- lapply(v, unlist)
        dir1 <- apply(dir1, 1, list)
        dir1 <- lapply(dir1, unlist)
        dir2 <- apply(dir2, 1, list)
        dir2 <- lapply(dir2, unlist)
        ans <- mapply(function(u, dir1){
          mapply(function(u, v, dir1, dir2){
            h <- A %*% matrix(u - v, 3, 1)
            d <- sqrt(sum(h^2))
            K <- exp(-3*d)
            h <- AtA %*% matrix(u - v, 3, 1)
            hess <- cont * K *
              (3*AtA/(d+eps) +
                 outer(h[,],h[,]) *
                 (9/(d^2+eps) +
                    3/(d^3+eps))
              )
            matrix(dir1, 1, 3) %*%
              hess %*%
              matrix(dir2, 3, 1)
          }, list(u), v, list(dir1), dir2)
        }, u, dir1)
        matrix(unlist(ans), length(u), length(v),
               byrow = T)
      }
      .Object@varfun <- function(u, v = u){
        .Object@contribution - .Object@covfun(u, v)
      }
      .Object@vard1 <- function(u, v, dir1){
        - .Object@covd1(u, v, dir1)
      }
      .Object@vard2 <- function(u, v, dir1, dir2){
        - .Object@covd2(u, v, dir1, dir2)
      }
    }
    else if (.Object@type == "spherical"){
      .Object@covfun <- function(u, v = u){
        A <- .Object@anis
        v1 <- t(A %*% t(v))
        u1 <- t(A %*% t(u))
        d <- .vectorized_pdist(u1, v1)
        d[d > 1] <- 1
        .Object@contribution * (1 - 1.5 * d + 0.5 * d^3)
      }
      .Object@covd1 <- function(u, v, dir1){
        eps <- .Machine$double.eps
        A <- .Object@anis
        AtA <- t(A) %*% A
        cont <- .Object@contribution
        u <- apply(u, 1, list)
        u <- lapply(u, unlist)
        v <- apply(v, 1, list)
        v <- lapply(v, unlist)
        dir1 <- apply(dir1, 1, list)
        dir1 <- lapply(dir1, unlist)
        ans <- lapply(u, function(u){
          mapply(function(u, v, dir1){
            h <- A %*% matrix(u - v, 3, 1)
            d <- min(1, sqrt(sum(h^2)))
            h <- matrix(u - v, 3, 1)
            grad <- cont * 1.5 *
              (1/(d+eps) - d) *
              AtA %*% h
            sum(grad * dir1)
          }, list(u), v, dir1)
        })
        matrix(unlist(ans), length(u), length(v),
               byrow = T)
      }
      .Object@covd2 <- function(u, v = u, dir1, dir2){
        eps <- .Machine$double.eps
        A <- .Object@anis
        AtA <- t(A) %*% A
        cont <- .Object@contribution
        u <- apply(u, 1, list)
        u <- lapply(u, unlist)
        v <- apply(v, 1, list)
        v <- lapply(v, unlist)
        dir1 <- apply(dir1, 1, list)
        dir1 <- lapply(dir1, unlist)
        dir2 <- apply(dir2, 1, list)
        dir2 <- lapply(dir2, unlist)
        ans <- mapply(function(u, dir1){
          mapply(function(u, v, dir1, dir2){
            h <- A %*% matrix(u - v, 3, 1)
            d <- min(1, sqrt(sum(h^2)))
            h <- AtA %*% matrix(u - v, 3, 1)
            hess <- cont * 1.5 * AtA * (1/(d + eps) - d) -
              1.5 * outer(h[,],h[,]) *
              (1/(d^3 + eps) + 1/(d + eps))
            matrix(dir1, 1, 3) %*%
              hess %*%
              matrix(dir2, 3, 1)
          }, list(u), v, list(dir1), dir2)
        }, u, dir1)
        matrix(unlist(ans), length(u), length(v),
               byrow = T)
      }
      .Object@varfun <- function(u, v = u){
        .Object@contribution - .Object@covfun(u, v)
      }
      .Object@vard1 <- function(u, v, dir1){
        - .Object@covd1(u, v, dir1)
      }
      .Object@vard2 <- function(u, v, dir1, dir2){
        - .Object@covd2(u, v, dir1, dir2)
      }
    }
    else if (.Object@type == "power"){
      # Obs.: the power model is hidden from the user in the documentation,
      # as it would require a special GP object that works with negative
      # definite variogram matrices. It may become available in a future
      # version.
      .Object@covfun <- function(u, v){
        .Object@contribution - .Object@varfun(u, v)
      }
      .Object@covd1 <- function(u, v, dir1){
        - .Object@vard1(u, v, dir1)
      }
      .Object@covd2 <- function(u, v, dir1, dir2){
        - .Object@vard2(u, v, dir1, dir2)
      }
      .Object@varfun <- function(u, v = u){
        A <- .Object@anis
        v1 <- t(A %*% t(v))
        u1 <- t(A %*% t(u))
        d <- .vectorized_pdist(u1, v1)
        .Object@contribution * d ^ (.Object@power)
      }
      .Object@vard1 <- function(u, v, dir1){
        eps <- .Machine$double.eps
        A <- .Object@anis
        AtA <- t(A) %*% A
        cont <- .Object@contribution
        p <- .Object@power
        u <- apply(u, 1, list)
        u <- lapply(u, unlist)
        v <- apply(v, 1, list)
        v <- lapply(v, unlist)
        dir1 <- apply(dir1, 1, list)
        dir1 <- lapply(dir1, unlist)
        ans <- lapply(u, function(u){
          mapply(function(u, v, dir1){
            h <- A %*% matrix(u - v, 3, 1)
            d <- sqrt(sum(h^2))
            h <- matrix(u - v, 3, 1)
            grad <- cont * p * d^p/(d^2+eps) *
              AtA %*% h
          }, list(u), v, dir1)
        })
        matrix(unlist(ans), length(u), length(v),
               byrow = T)
      }
      .Object@vard2 <- function(u, v = u, dir1, dir2){
        eps <- .Machine$double.eps
        A <- .Object@anis
        AtA <- t(A) %*% A
        cont <- .Object@contribution
        u <- apply(u, 1, list)
        u <- lapply(u, unlist)
        v <- apply(v, 1, list)
        v <- lapply(v, unlist)
        dir1 <- apply(dir1, 1, list)
        dir1 <- lapply(dir1, unlist)
        dir2 <- apply(dir2, 1, list)
        dir2 <- lapply(dir2, unlist)
        ans <- mapply(function(u, dir1){
          mapply(function(u, v, dir1, dir2){
            h <- A %*% matrix(u - v, 3, 1)
            d <- sqrt(sum(h^2))
            h <- AtA %*% matrix(u - v, 3, 1)
            hess <- cont * (p*(p-2) *
                              d^p/(d^4+eps) *
                              outer(h[,],h[,]) +
                              p * d^p/(d^2+eps) * AtA)
            matrix(dir1, 1, 3) %*%
              hess %*%
              matrix(dir2, 3, 1)
          }, list(u), v, list(dir1), dir2)
        }, u, dir1)
        matrix(unlist(ans), length(u), length(v),
               byrow = T)
      }
    }
    else if (.Object@type == "cubic"){
      .Object@covfun <- function(u, v = u){
        A <- .Object@anis
        v1 <- t(A %*% t(v))
        u1 <- t(A %*% t(u))
        d <- .vectorized_pdist(u1, v1)
        d[d > 1] <- 1
        .Object@contribution * (1 - 7 * d^2 + 35/4 * d^3 - 7/2 * d^5 +
                                  3/4 * d^7)
      }
      .Object@covd1 <- function(u, v, dir1){
        .Object@contribution *
          .covd1_cubic(u, v, dir1, .Object@anis)
      }
      .Object@covd2 <- function(u, v = u, dir1, dir2){
        .Object@contribution *
          .covd2_cubic(u, v, dir1, dir2,
                      .Object@anis)
      }
      .Object@varfun <- function(u, v = u){
        .Object@contribution - .Object@covfun(u, v)
      }
      .Object@vard1 <- function(u, v, dir1){
        - .Object@covd1(u, v, dir1)
      }
      .Object@vard2 <- function(u, v, dir1, dir2){
        - .Object@covd2(u, v, dir1, dir2)
      }
    }
    else if (.Object@type == "matern1"){
      .Object@covfun <- function(u, v = u){
        A <- .Object@anis
        v1 <- t(A %*% t(v))
        u1 <- t(A %*% t(u))
        d <- .vectorized_pdist(u1, v1)
        .Object@contribution * (1 + 5 * d) * exp(-5 * d)
      }
      .Object@covd1 <- function(u, v, dir1){
        .Object@contribution *
          .covd1_matern1(u, v, dir1, .Object@anis)
      }
      .Object@covd2 <- function(u, v = u, dir1, dir2){
        .Object@contribution *
          .covd2_matern1(u, v, dir1, dir2, .Object@anis)
      }
      .Object@varfun <- function(u, v = u){
        .Object@contribution - .Object@covfun(u, v)
      }
      .Object@vard1 <- function(u, v, dir1){
        - .Object@covd1(u, v, dir1)
      }
      .Object@vard2 <- function(u, v, dir1, dir2){
        - .Object@covd2(u, v, dir1, dir2)
      }
    }
    else if (.Object@type == "matern2"){
      .Object@covfun <- function(u, v = u){
        A <- .Object@anis
        v1 <- t(A %*% t(v))
        u1 <- t(A %*% t(u))
        d <- .vectorized_pdist(u1, v1)
        .Object@contribution * (1 + 5*d + 12*d^2) * exp(-6 * d)
      }
      .Object@covd1 <- function(u, v, dir1){
        .Object@contribution *
          .covd1_matern2(u, v, dir1, .Object@anis)
      }
      .Object@covd2 <- function(u, v = u, dir1, dir2){
        .Object@contribution *
          .covd2_matern2(u, v, dir1, dir2, .Object@anis)
      }
      .Object@varfun <- function(u, v = u){
        .Object@contribution - .Object@covfun(u, v)
      }
      .Object@vard1 <- function(u, v, dir1){
        - .Object@covd1(u, v, dir1)
      }
      .Object@vard2 <- function(u, v, dir1, dir2){
        - .Object@covd2(u, v, dir1, dir2)
      }
    }
    else if (.Object@type == "cauchy"){
      .Object@covfun <- function(u, v){
        A <- .Object@anis
        v1 <- t(A %*% t(v))
        u1 <- t(A %*% t(u))
        beta <- 0.05^(-1/.Object@power) - 1
        d <- .vectorized_pdist(u1, v1)
        .Object@contribution * (1+beta*d^2) ^ (- .Object@power)
      }
      .Object@covd1 <- function(u, v, dir1){
        .Object@contribution *
          .covd1_cauchy(u, v, dir1, .Object@anis, .Object@power)
      }
      .Object@covd2 <- function(u, v, dir1, dir2){
        .Object@contribution *
          .covd2_cauchy(u, v, dir1, dir2, .Object@anis, .Object@power)
      }
      .Object@varfun <- function(u, v = u){
        .Object@contribution - .Object@covfun(u, v)
      }
      .Object@vard1 <- function(u, v, dir1){
        - .Object@covd1(u, v, dir1)
      }
      .Object@vard2 <- function(u, v, dir1, dir2){
        - .Object@covd2(u, v, dir1, dir2)
      }
    }
    return(.Object)
  }
)

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

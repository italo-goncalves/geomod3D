#' @include generics.R
#' @include utils.R
NULL

#### covarianceStructure3D class ####
#' Covariance structure
#'
#' Representation of one structure in a spatial continuity model.
#'
#' @slot type The type of covariance function represented.
#' @slot params A list containing the specific parameters for each covariance
#' type.
#'
#' @seealso \code{\link{covarianceStructure3D-init}},
#' \code{\link{CovarianceMatrix}}
#'
#' @name covarianceStructure3D-class
#' @export covarianceStructure3D
covarianceStructure3D <- setClass(
  "covarianceStructure3D",
  slots = c(type = "character",
            params = "list"),
  validity = function(object) {
    if (!(object@type %in% c("gaussian",
                             "exponential",
                             "spherical",
                             "cubic",
                             "matern1",
                             "matern2",
                             "cauchy",
                             "linear",
                             "bias",
                             "nnet")))
      stop(paste0("Structure type '",object@type, "' not supported."))
    if (object@params$contribution < 0)
      stop("Contribution cannot be negative")
    if (min(c(object@params$maxrange, object@params$midrange,
              object@params$minrange)) <= 0)
      stop("All ranges must be greater than zero")
    if (object@type == "cauchy" && object@params$power <= 0)
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
#' @param center A vector with the XYZ coordinates for centering some
#' covariance functions (see details).
#'
#' @section Covariance functions:
#'
#' The following functions are available for the \code{type} parameter. In the
#' equations below, \eqn{X1} and \eqn{X2} are coordinate vectors,
#' \eqn{d} is the anisotropy-corrected distance between two points,
#' \eqn{p} is the \code{power} parameter, and \eqn{c} is the center. Functions
#' marked with an asterisk support modeling with structural data.
#'
#' \describe{
#'   \item{\code{gaussian}*}{\deqn{\exp{-3d^2}}{exp(-3*d^2)}}
#'   \item{\code{exponential}}{\deqn{\exp{-3d}}{exp(-3*d)}}
#'   \item{\code{spherical}}{\deqn{1 - 1.5 d + 0.5 d^3, 0 \leq d \leq 1}{1 - 1.5*d + 0.5*d^3, 0 <= d <= 1}}
#'   \item{\code{cubic}*}{\deqn{1 - 7d^2 + \frac{35}{4}d^3 - \frac{7}{2}d^5 + \frac{3}{4}d^7, 0 \leq d \leq 1}{1 - 7*d^2 + 35/4*d^3 - 7/2*d^5 + 3/4*d^7, 0 <= d <= 1}}
#'   \item{\code{matern1}*}{\deqn{(1 + 5d)\exp{-5d}}{(1 + 5*d) * exp(-5*d)}}
#'   \item{\code{matern2}*}{\deqn{(1 + 6d + 12d^2)\exp{-6d}}{(1 + 6*d + 12*d^2) * exp(-6*d)}}
#'   \item{\code{cauchy}*}{\deqn{(1 + d^2)^{-p}, p > 0}{(1 + d^2)^(-p), p > 0}}
#'   \item{\code{linear}*}{\deqn{(X1 - c)^T(X2 - c)}{(X1 - c)^T(X2 - c)}}
#'   \item{\code{bias}}{\deqn{1}{1}}
#' }
#'
#' @name covarianceStructure3D-init
#'
#' @seealso \code{\link{covarianceStructure3D-class}},
#' \code{\link{covarianceModel3D-class}}
covarianceStructure3D <- function(type = c("gaussian", "exponential",
                                           "spherical", "cubic", "matern1",
                                           "matern2", "cauchy",
                                           "linear", "bias", "nnet"),
                                  contribution = 1,
                                  maxrange = 1, midrange = maxrange,
                                  minrange = midrange,
                                  azimuth = 0, dip = 0, rake = 0,
                                  power = 1, center = c(0, 0, 0)){
  type <- type[1]

  params <- list(contribution = contribution)

  if (type != "bias")
  params <- c(params, list(maxrange = maxrange,
                           midrange = midrange,
                           minrange = minrange,
                           azimuth = azimuth,
                           dip = dip,
                           rake = rake))

  if (type %in% c("cauchy", "nnet"))
    params <- c(params, list(power = power))

  if (type %in% c("linear", "nnet"))
    params <- c(params, list(center = center))

  new("covarianceStructure3D", type = type, params = params)
}

#### show ####
setMethod(
  f = "show",
  signature = "covarianceStructure3D",
  definition = function(object){
    cat("Object of class ", class(object), "\n\n", sep = "")
    cat("Type:", object@type,"\n")
    cat("Contribution =", object@params$contribution,"\n")
    if (type != "bias"){
      cat("Maximum range =", object@params$maxrange,"\n")
      cat("Intermediate range =", object@params$midrange,"\n")
      cat("Minimum range =", object@params$minrange,"\n")
      cat("Orientation: azimuth =", object@params$azimuth,
          "dip =", object@params$dip, "rake =", object@params$rake, "\n")
    }

    if(object@type %in% c("cauchy", "nnet")){
      cat("Power =", object@params$power,"\n")
    }
    if (object@type %in% c("linear", "nnet")) {
      center <- object@params$center
      cat("Center (X Y Z) = ", center,"\n")
    }
  }
)

#### Covariance Matrix ####
setMethod(
  f = ".CalcCovMat",
  signature = "covarianceStructure3D",
  definition = function(object, u, v){
    if (object@type == "bias")
      return(matrix(object@params$contribution), nrow(u), nrow(v))

    if (object@type %in% c("linear", "nnet")){
      u <- u - matrix(object@params$center, nrow(u), 3, byrow = T)
      v <- v - matrix(object@params$center, nrow(v), 3, byrow = T)
    }

    A <- solve(anisotropy3D(object@params$maxrange,
                            object@params$midrange,
                            object@params$minrange,
                            object@params$azimuth,
                            object@params$dip,
                            object@params$rake))
    contribution <- object@params$contribution
    v1 <- t(A %*% t(v))
    u1 <- t(A %*% t(u))

    if (object@type == "gaussian"){
      d <- .vectorized_pdist(u1, v1)
      contribution * exp(-3*(d^2))
    }
    else if (object@type == "exponential"){
      d <- .vectorized_pdist(u1, v1)
      contribution * exp(-3 * d)
    }
    else if (object@type == "spherical"){
      d <- .vectorized_pdist(u1, v1)
      d[d > 1] <- 1
      contribution * (1 - 1.5 * d + 0.5 * d^3)
    }
    else if (object@type == "cubic"){
      d <- .vectorized_pdist(u1, v1)
      d[d > 1] <- 1
      contribution * (1 - 7 * d^2 + 35/4 * d^3 - 7/2 * d^5 + 3/4 * d^7)
    }
    else if (object@type == "matern1"){
      d <- .vectorized_pdist(u1, v1)
      contribution * (1 + 5 * d) * exp(-5 * d)
    }
    else if (object@type == "matern2"){
      d <- .vectorized_pdist(u1, v1)
      contribution * (1 + 5*d + 12*d^2) * exp(-6 * d)
    }
    else if (object@type == "cauchy"){
      d <- .vectorized_pdist(u1, v1)
      p <- object@params$power
      beta <- 0.05^(-1/p) - 1
      contribution * (1+beta*d^2) ^ (- p)
    }
    else if (object@type == "linear"){
      contribution * tcrossprod(u, v)
    }
    else if (object@type == "nnet"){
      A <- solve(A)

      sigma <- Diagonal(4, 1);
      sigma[1, 1] <- object@params$power
      sigma[2:4, 2:4] <- A
      u <- cbind(1, u); v <- cbind(1, v)

      tmp <- crossprod(t(u), sigma %*% t(v)) * 2
      tmp1 <- colSums(t(u) * (sigma %*% t(u))) * 2 + 1
      tmp1 <- Matrix(tmp1, nrow(tmp), ncol(tmp))
      tmp2 <- colSums(t(v) * (sigma %*% t(v))) * 2 + 1
      tmp2 <- Matrix(tmp2, nrow(tmp), ncol(tmp), byrow = T)

      tmp3 <- sqrt(tmp1 * tmp2)
      K <- asin(tmp / tmp3)

      K * contribution * 2 / pi
    }
  }
)

setMethod(
  f = ".CalcCovMat_d1",
  signature = "covarianceStructure3D",
  definition = function(object, u, v, dir1){
    if (object@type %in% c("exponential", "spherical", "bias"))
      stop(paste(object@type, "function does not support directional data"))

    anis <- solve(anisotropy3D(object@params$maxrange,
                               object@params$midrange,
                               object@params$minrange,
                               object@params$azimuth,
                               object@params$dip,
                               object@params$rake))
    contribution <- object@params$contribution

    if (object@type == "gaussian"){
      contribution * .covd1_gaussian(u, v, dir1, anis)
    }
    else if (object@type == "cubic"){
      contribution * .covd1_cubic(u, v, dir1, anis)
    }
    else if (object@type == "matern1"){
      contribution * .covd1_matern1(u, v, dir1, anis)
    }
    else if (object@type == "matern2"){
      contribution * .covd1_matern2(u, v, dir1, anis)
    }
    else if (object@type == "cauchy"){
      contribution * .covd1_cauchy(u, v, dir1, anis, object@params$power)
    }
    else if (object@type == "linear"){
      A <- crossprod(anis)
      contribution * tcrossprod(u %*% A, dir1)
    }
    else if (object@type == "nnet"){
      stop("nnet not yet implemented")
    }
  }
)

setMethod(
  f = ".CalcCovMat_d2",
  signature = "covarianceStructure3D",
  definition = function(object, u, v, dir1, dir2){
    if (object@type %in% c("exponential", "spherical", "bias"))
      stop(paste(object@type, "function does not support directional data"))

    anis <- solve(anisotropy3D(object@params$maxrange,
                               object@params$midrange,
                               object@params$minrange,
                               object@params$azimuth,
                               object@params$dip,
                               object@params$rake))
    contribution <- object@params$contribution

    if (object@type == "gaussian"){
      contribution * .covd2_gaussian(u, v, dir1, dir2, anis)
    }
    else if (object@type == "exponential"){
      stop("exponential function does not support directional data")
    }
    else if (object@type == "spherical"){
      stop("spherical function does not support directional data")
    }
    else if (object@type == "cubic"){
      contribution * .covd2_cubic(u, v, dir1, dir2, anis)
    }
    else if (object@type == "matern1"){
      contribution * .covd2_matern1(u, v, dir1, dir2, anis)
    }
    else if (object@type == "matern2"){
      contribution * .covd2_matern2(u, v, dir1, dir2, anis)
    }
    else if (object@type == "cauchy"){
      contribution * .covd2_cauchy(u, v, dir1, dir2, anis, object@power)
    }
    else if (object@type == "linear"){
      A <- crossprod(anis)
      contribution * tcrossprod(dir1 %*% A, dir2)
    }
    else if (object@type == "nnet"){
      stop("nnet not yet implemented")
    }
  }
)

#' @rdname GetPriorVariance
setMethod(
  f = "GetPriorVariance",
  signature = "covarianceStructure3D",
  definition = function(object, target){
    if (object@type == "linear"){
      coords <- GetCoords(target, "matrix")
      coords <- coords - matrix(object@params$center, nrow(target), 3, byrow = T)
      return(object@params$contribution * rowSums(coords * coords))
    }
    else
      return(rep(object@params$contribution, nrow(target)))
  }
)

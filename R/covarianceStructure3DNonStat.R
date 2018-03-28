#' @include generics.R
#' @include utils.R
NULL

#### covarianceStructure3D class ####
#' Non stationary covariance structure
#'
#' Representation of one structure in a non stationary spatial continuity model.
#'
#' @slot type The type of covariance function represented.
#' @slot contribution The variance associated with this structure.
#' @slot maxrange,midrange,minrange The structure's range in three orthogonal
#' directions, or the three semi-axes of the covariance ellipsoid.
#' @slot azimuth,dip,rake Orientation of the covariance ellipsoid in geological
#' coordinates.
#' @slot power The power parameter for the covariance functions that feature
#' one.
#' @slot covfun,covd1,covd2 Functions to calculate the
#' covariance matrices between a variable of interest and/or its
#' derivatives.
#'
#' @seealso \code{\link{covarianceStructure3DNonStat-init}},
#' \code{\link{CovarianceMatrix}}
#'
#' @name covarianceStructure3DNonStat-class
#' @export covarianceStructure3DNonStat
covarianceStructure3DNonStat <- setClass(
  "covarianceStructure3DNonStat",
  slots = c(type = "character",
            contribution = "GP",
            maxrange = "GP",
            midrange = "GP",
            minrange = "GP",
            azimuth = "GP",
            dip = "GP",
            rake = "GP",
            power = "numeric"),
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
    if (object@type == "cauchy" && object@power <= 0)
      stop("Power must be positive for cauchy function")
    return(TRUE)
  }
)

#### initialize ####
#' Non stationary covariance structure
#'
#' Representation of one structure in a non stationary spatial continuity model.
#'
#' @param type The type of covariance function represented (see below).
#' @param contribution The variance associated with this structure.
#' @param maxrange,midrange,minrange The structure's range in three orthogonal
#' directions, or the three semi-axes of the covariance ellipsoid.
#' @param azimuth,dip,rake Orientation of the covariance ellipsoid in geological
#' coordinates.
#' @param power The power parameter for the covariance functions that feature
#' one.
#' @param control_points A \code{points3DDataFrame} object with the coordinates
#' of the non stationary features.
#'
#' @details All parameters except \code{power} are non stationary, each
#' modeled by its own (stationary, noiseless) \code{GP} and using the control points as
#' data. The GPs are initialized with the values provided as data, and by
#' default predict a constant value everywhere. Use the \code{Fit} method of the
#' main GP to determine the best value for each parameter at the control points
#' and the stationary covariance model for each.
#'
#' @section Covariance functions:
#'
#' The following functions are available for the \code{type} parameter. In the
#' equations below, \eqn{d} is the anisotropy-corrected distance between two points
#' and \eqn{p} is the \code{power} parameter. Functions marked with an asterisk
#' support modeling with directional data.
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
#' @name covarianceStructure3DNonStat-init
#'
#' @seealso \code{\link{covarianceStructure3DNonStat-class}},
#' \code{\link{covarianceModel3D-class}}
covarianceStructure3DNonStat <- function(type = c("gaussian", "exponential",
                                           "spherical", "cubic", "matern1",
                                           "matern2", "cauchy"),
                                  contribution,
                                  maxrange, midrange = maxrange,
                                  minrange = midrange,
                                  azimuth = 0, dip = 0, rake = 0,
                                  power = 1,
                                  control_points,
                                  pseudo_inputs = control_points){
  type <- type[1]

  # contribution and maxrange are modeled as log to ensure positivity
  control_points[, "contribution"] <- log(contribution)
  control_points[, "maxrange"] <- log(maxrange)

  # midrange and minrange are modeled as logit to bound between 0.05 and 1
  # (midrange is modeled as a fraction of maxrange, and
  # minrange is modeled as a fraction of midrange)
  control_points[, "midrange"] <- .logit(min(0.999, midrange / maxrange))
  control_points[, "minrange"] <- .logit(min(0.999, minrange / midrange))

  # change of angular system
  if (azimuth >= 180){
    azimuth <- azimuth - 180
    dip <- -dip
  }

  # angles also modeled as logit
  control_points[, "azimuth"] <- .logit(min(0.999, max(0.001, azimuth / 180)))
  control_points[, "dip"] <- .logit(min(0.999, max(0.001, (dip + 90) / 180)))
  control_points[, "rake"] <- .logit(min(0.999, max(0.001, (rake + 90) / 180)))

  # initializing GPs
  contrGP <- SPGP(control_points, value = "contribution", reg.v = 1e-5,
                model = covarianceModel3D(0,
                          covarianceStructure3D("gaussian", 0, 1)),
                pseudo_inputs = pseudo_inputs)
  maxrGP <- SPGP(control_points, value = "maxrange", reg.v = 1e-5,
                model = covarianceModel3D(0,
                                covarianceStructure3D("gaussian", 0, 1)),
                pseudo_inputs = pseudo_inputs)
  midrGP <- SPGP(control_points, value = "midrange", reg.v = 1e-5,
               model = covarianceModel3D(0,
                                covarianceStructure3D("gaussian", 0, 1)),
               pseudo_inputs = pseudo_inputs)
  minrGP <- SPGP(control_points, value = "minrange", reg.v = 1e-5,
               model = covarianceModel3D(0,
                              covarianceStructure3D("gaussian", 0, 1)),
               pseudo_inputs = pseudo_inputs)
  azGP <- SPGP(control_points, value = "azimuth", reg.v = 1e-5,
               model = covarianceModel3D(0,
                            covarianceStructure3D("gaussian", 0, 1)),
               pseudo_inputs = pseudo_inputs)
  dipGP <- SPGP(control_points, value = "dip", reg.v = 1e-5,
               model = covarianceModel3D(0,
                            covarianceStructure3D("gaussian", 0, 1)),
               pseudo_inputs = pseudo_inputs)
  rakeGP <- SPGP(control_points, value = "rake", reg.v = 1e-5,
               model = covarianceModel3D(0,
                          covarianceStructure3D("gaussian", 0, 1)),
               pseudo_inputs = pseudo_inputs)

  # end
  new("covarianceStructure3DNonStat", type = type, contribution = contrGP,
      maxrange = maxrGP, midrange = midrGP, minrange = minrGP,
      azimuth = azGP, dip = dipGP, rake = rakeGP, power = power)
}

#### show ####
setMethod(
  f = "show",
  signature = "covarianceStructure3DNonStat",
  definition = function(object){
    cat("Object of class ", class(object), "\n\n", sep = "")
    cat("Type:", object@type,"\n")
  }
)

#### logLik ####
setMethod(
  f = "logLik",
  signature = "covarianceStructure3DNonStat",
  definition = function(object){
    sum(
      logLik(object@contribution),
      logLik(object@maxrange),
      logLik(object@midrange),
      logLik(object@minrange),
      logLik(object@azimuth),
      logLik(object@dip),
      logLik(object@rake)
    )
  }
)

#### GetNonStatParams ####
#' @rdname GetNonStatParams
setMethod(
  f = "GetNonStatParams",
  signature = "covarianceStructure3DNonStat",
  definition = function(object, target){
    target <- Predict(object@contribution, target, to = "contribution", output.var = F)
    target[, "contribution"] <- exp(target[["contribution"]])

    target <- Predict(object@maxrange, target, to = "maxrange", output.var = F)
    target[, "maxrange"] <- exp(target[["maxrange"]])

    target <- Predict(object@midrange, target, to = "midrange", output.var = F)
    target[, "midrange"] <- (0.05 + 0.95 * .sig(target[["midrange"]])) *
      target[["maxrange"]]

    target <- Predict(object@minrange, target, to = "minrange", output.var = F)
    target[, "minrange"] <- (0.05 + 0.95 * .sig(target[["minrange"]])) *
      target[["midrange"]]

    target <- Predict(object@azimuth, target, to = "azimuth", output.var = F)
    target[, "azimuth"] <- .sig(target[["azimuth"]]) * 180

    target <- Predict(object@dip, target, to = "dip", output.var = F)
    target[, "dip"] <- .sig(target[["dip"]]) * 180 - 90

    target <- Predict(object@rake, target, to = "rake", output.var = F)
    target[, "rake"] <- .sig(target[["rake"]]) * 180 - 90

    flip <- target[["dip"]] < 0
    target[flip, "dip"] <- - target[["dip"]][flip]
    target[flip, "azimuth"] <- target[["azimuth"]][flip] + 180

    return(target)
  }
)


#### Covariance Matrix ####
setMethod(
  f = ".CalcCovMat",
  signature = "covarianceStructure3DNonStat",
  definition = function(object, u, v){

    # interpolation and conversion
    u3d <- points3DDataFrame(u, data.frame(contribution = rep(NA, nrow(u))))
    u3d <- GetNonStatParams(object, u3d)

    v3d <- points3DDataFrame(v, data.frame(contribution = rep(NA, nrow(v))))
    v3d <- GetNonStatParams(object, v3d)

    # output
    K <- .cov_ns(u, v, sqrt(u3d[["contribution"]]), sqrt(v3d[["contribution"]]),
            u3d[["maxrange"]], v3d[["maxrange"]],
            u3d[["midrange"]], v3d[["midrange"]],
            u3d[["minrange"]], v3d[["minrange"]],
            u3d[["azimuth"]], v3d[["azimuth"]],
            u3d[["dip"]], v3d[["dip"]],
            u3d[["rake"]], v3d[["rake"]],
            object@type, object@power)
    if (any(is.nan(K))) browser()

    return(K)

  }
)

# setMethod(
#   f = ".CalcCovMat_d1",
#   signature = "covarianceStructure3D",
#   definition = function(object, u, v, dir1){
#     anis <- object@anis
#     contribution <- object@contribution
#
#     if (object@type == "gaussian"){
#       contribution * .covd1_gaussian(u, v, dir1, anis)
#     }
#     else if (object@type == "exponential"){
#       stop("exponential function does not support directional data")
#     }
#     else if (object@type == "spherical"){
#       stop("spherical function does not support directional data")
#     }
#     else if (object@type == "cubic"){
#       contribution * .covd1_cubic(u, v, dir1, anis)
#     }
#     else if (object@type == "matern1"){
#       contribution * .covd1_matern1(u, v, dir1, anis)
#     }
#     else if (object@type == "matern2"){
#       contribution * .covd1_matern2(u, v, dir1, anis)
#     }
#     else if (object@type == "cauchy"){
#       contribution * .covd1_cauchy(u, v, dir1, anis, object@power)
#     }
#   }
# )
#
# setMethod(
#   f = ".CalcCovMat_d2",
#   signature = "covarianceStructure3D",
#   definition = function(object, u, v, dir1, dir2){
#     anis <- object@anis
#     contribution <- object@contribution
#
#     if (object@type == "gaussian"){
#       contribution * .covd2_gaussian(u, v, dir1, dir2, anis)
#     }
#     else if (object@type == "exponential"){
#       stop("exponential function does not support directional data")
#     }
#     else if (object@type == "spherical"){
#       stop("spherical function does not support directional data")
#     }
#     else if (object@type == "cubic"){
#       contribution * .covd2_cubic(u, v, dir1, dir2, anis)
#     }
#     else if (object@type == "matern1"){
#       contribution * .covd2_matern1(u, v, dir1, dir2, anis)
#     }
#     else if (object@type == "matern2"){
#       contribution * .covd2_matern2(u, v, dir1, dir2, anis)
#     }
#     else if (object@type == "cauchy"){
#       contribution * .covd2_cauchy(u, v, dir1, dir2, anis, object@power)
#     }
#   }
# )

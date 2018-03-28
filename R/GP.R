#' @include generics.R
#' @include utils.R
#' @include spatial3DDataFrame.R
#' @include points3DDataFrame.R
#' @include directions3DDataFrame.R
#' @include covarianceModel3D.R
NULL

#### Gaussian Process class ####
#' Gaussian Process
#'
#' Implementation of the Gaussian Process model for 3D spatial interpolation.
#'
#' @slot data A \code{spatial3DDataFrame} object containing the necessary data.
#' @slot tangents A \code{directions3DDataFrame} object containing structural
#' geology data. Most likely generated with the \code{GetLineDirections()}
#' or \code{GetPlaneDirections()} method.
#' @slot model The covariance model. A \code{list} containing one or more
#' \code{covarianceStructure3D} objects.
#' @slot mean The global mean. Irrelevant if a trend is used.
#' @slot trend The model's trend component. A formula in character format.
#' @slot beta The trend coefficients.
#' @slot likelihood The model's log-likelihood given the data.
#' @slot pre_comp A \code{list} containing pre-computed values to speed up
#' predictions.
#'
#' @references Rasmussen CE, Williams CKI. Gaussian processes for machine
#' learning. Cambridge, Massachusetts: MIT Press; 2006.
#' doi:10.1142/S0129065704001899.
#'
#' @seealso \code{\link{GP-init}}, \code{\link{SPGP-class}},
#' \code{\link{GP_geomod-class}}
#'
#' @name GP-class
#' @export GP
GP <- setClass(
  "GP",
  slots = c(data = "spatial3DDataFrame",
            tangents = "directions3DDataFrame",
            model = "covarianceModel3D",
            mean = "numeric",
            trend = "character",
            beta = "matrix",
            likelihood = "numeric",
            pre_comp = "list")
)

#### initialization ####
#' Gaussian Process
#'
#' Implementation of the Gaussian Process model for 3D spatial interpolation.
#'
#' @param data A \code{spatial3DDataFrame} object containing the data one
#' wishes to model.
#' @param model The covariance model. A \code{covarianceModel3D} object.
#' @param value The column name of the variable to be modeled. It is assumed
#' the column does not contain missing values.
#' @param mean The global mean. Irrelevant if a trend is provided.
#' @param trend The model's trend component. A formula in character format.
#' @param weights The importance of each data point in the model (a vector with
#' values between 0 and 1)
#' @param force.interp Indices of points that must be interpolated exactly.
#' @param reg.v Regularization to improve stability. A single value or a vector
#' with length matching the number of data points.
#' @param tangents A \code{directions3DDataFrame} object containing structural
#' geology data. Most likely generated with the \code{GetLineDirections()}
#' method.
#' @param reg.t Regularization for structural data. A single value or a vector
#' with length matching the number of structural data.
#'
#' @details This method builds a \code{GP} object with all the information
#' needed to make preditions at new data points.
#'
#' \code{trend} must be a character string with a formula as a function of
#' uppercase X, Y, and Z. The most common is the linear trend,
#' \code{"~ X + Y + Z"}. For ordinary kriging, use \code{"~1"}. If neither
#' \code{trend} nor \code{mean} are given, it is assumed that the global mean
#' is the mean of the data values.
#'
#' If any point index is given in \code{force.interp}, the predicted mean
#' function will pass exactly through those points, but the predictive variance
#' will still be computed as usual. This is in contrast to what is usually
#' done by geostatistics softwares, which assign a variance of 0 to those
#' points.
#'
#' \code{weights} can be used to "tune down" the effect of some data points.
#' The smaller the weight, the less effect a point will have on the predicted
#' funtion, to the limit that a weight of 0 filters the point completely.
#'
#' Note that this implementation uses all the data provided to make predictions,
#' which may be too memory intensive for large datasets.
#'
#' @name GP-init
#'
#' @seealso \code{\link{GP-class}}, \code{\link{SPGP-class}}
GP <- function(data, model, value,
               mean = NULL, trend = NULL,
               force.interp = numeric(), reg.v = 1e-9,
               tangents = NULL, reg.t = 1e-12, nugget.t = 0){

  # value
  if(length(value) == 1 & class(value) == "character")
    data["value"] <- data[value]
  else
    data["value"] <- value
  yval <- data[["value"]]

  # interpolation
  int <- rep(F, nrow(data))
  int[force.interp] <- T
  data["interpolate"] <- int

  # regularization for tangents
  # if(!is.null(tangents) && nrow(tangents) > 0){
  # if(length(reg.t) == 1)
  # tangents["reg"] <- rep(reg.t, nrow(tangents))
  #   else
  #     tangents["reg"] <- reg.t
  #   reg.t <- tangents[["reg"]]
  # }

  # mean
  if(is.null(mean)) mean <- mean(yval)

  # data
  data2 <- data[c("value", "interpolate")]
  tangents <- as(tangents, "directions3DDataFrame")

  # trend
  if(is.null(trend) | length(trend) == 0){
    TR <- matrix(0, nrow(data), 0)
    trend <- character(0)
  }else{
    TR <- TrendMatrix(data, trend)
    if(!is.null(tangents) && nrow(tangents) > 0){
      TR <- rbind(
        TR,
        TrendMatrix(tangents, trend)
      )
    }
    mean <- 0
  }
  Ntrend <- dim(TR)[2]

  # covariances
  Ntang <- 0
  K <- CovarianceMatrix(data, data, model)
  # nugget
  nug2 <-  rep(model@nugget, nrow(data))
  nug2[force.interp] <- reg.v
  K <- K + diag(nug2, length(nug2), length(nug2))
  if(!is.null(tangents) && nrow(tangents) > 0){
    K1 <- CovarianceMatrix(data, tangents, model)
    K2 <- CovarianceMatrix(tangents, tangents, model)
    Ntang <- nrow(K2)
    # regularization
    K2 <- K2 + diag(reg.t + model@nugget.dir, Ntang, Ntang)
    # final matrix
    K <- rbind(
      cbind(K, K1),
      cbind(t(K1), K2)
    )
  }
  # enforcing symmetry
  # K <- 0.5 * K + 0.5 * t(K)

  # pre-computations
  pre_comp <- list()
  yval <- c(yval - mean, rep(0, Ntang))

  # L <- t(chol(Matrix(K)))  # makes L lower triangular
  L <- .safeChol(Matrix(K))
  LiY <- solve(L, Matrix(yval, length(yval), 1))
  w_value <- solve(t(L), LiY)

  pre_comp$w_value <- as.numeric(w_value)
  pre_comp$w_var <- L
  beta <- matrix(0, 0, 0)
  if(Ntrend > 0){
    HLi <- t(TR) %*% solve(L)
    A <- HLi %*% t(HLi)
    b1 <- t(TR) %*% w_value
    beta <- solve(A, b1)
    beta <- as.matrix(beta)
    pre_comp$w_trend <- HLi
    pre_comp$A <- A
  }

  # likelihood
  likelihood <- - sum(log(diag(L))) - 0.5 * sum(yval * w_value) -
    0.5 * length(yval) * log(2 * pi)
  if(Ntrend > 0){
    LiYH <- HLi %*% LiY
    tmp <- t(LiYH) %*% solve(A, LiYH)
    likelihood <- likelihood + 0.5 * as.numeric(tmp) +
      Ntrend * log(2*pi) - 0.5 * determinant(A)$modulus
  }

  # end
  new("GP", data = data2, tangents = tangents, model = model, mean = mean,
      trend = trend, beta = beta, likelihood = likelihood, pre_comp = pre_comp)
}

#### show ####
setMethod(
  f = "show",
  signature = "GP",
  definition = function(object){
    # display
    cat("Object of class ", class(object), "\n", sep = "")
    cat("Data points:", nrow(object@data), "\n")
    cat("Tangent points:", nrow(object@tangents), "\n")
    if(length(object@trend) == 0)
      cat("Global mean:", object@mean, "\n")
    else{
      cat("Trend:", object@trend, "\n")
      print(object@beta)
    }
    cat("Log-likelihood:", object@likelihood, "\n")
  }
)

#### summary ####
setMethod(
  f = "summary",
  signature = "GP",
  definition = function(object){
    show(object)
    cat("\n")
    show(object@model)

    # tangents
    if (nrow(object@tangents) > 0)
      cat("Nugget for tangents:", object@model@nugget.t, "\n")
  }
)


#### Predict ####
#' @rdname Predict
setMethod(
  f = "Predict",
  signature = "GP",
  definition = function(object, target, to = "value", output.var = T){

    prior.var <- GetPriorVariance(object@model, target)

    # pre processing
    w_var <- Matrix(object@pre_comp$w_var)
    w_value <- object@pre_comp$w_value

    # trend
    if(length(object@trend) > 0){
      TRtarget <- TrendMatrix(target, object@trend)
      TRdata <- TrendMatrix(object@data, object@trend)
      w_tr <- object@pre_comp$w_trend
      beta <- object@beta
    }

    # trivial solution
    if(all(prior.var == 0)){
      target[, to] <- object@mean
      if(output.var) target[, paste0(to, ".var")] <- 0
      return(target)
    }

    # covariances
    Ntang <- nrow(object@tangents)
    Ktarget <- Matrix(CovarianceMatrix(target, object@data, object@model))
    if(Ntang > 0){
      K1 <- Matrix(CovarianceMatrix(target, object@tangents, object@model))
      Ktarget <- cbind(Ktarget, K1)
    }

    # prediction
    # residuals
    pred <- apply(Ktarget, 1, function(rw) sum(rw * w_value)) + object@mean
    # trend
    if(length(object@trend) > 0){
      LinvK <- solve(w_var, t(Ktarget))
      R <- t(TRtarget) - w_tr %*% LinvK
      pred <- pred + t(R) %*% beta
    }
    target[, to] <- pred

    # variance
    if(output.var){
      # tot_var <- object@model@total.var
      if(length(object@trend) > 0){
        pred_var <- colSums(LinvK ^ 2)
        pred_var[pred_var > prior.var] <- prior.var
        pred_var <- prior.var - pred_var + object@model@nugget
        tr_var <- colSums(
          R * (solve(w_tr %*% t(w_tr), R))
        )
        pred_var <- pred_var + tr_var
      }
      else{
        pred_var <- colSums(solve(w_var, t(Ktarget)) ^ 2)
        pred_var[pred_var > prior.var] <- prior.var
        pred_var <- prior.var - pred_var + object@model@nugget
      }
      target[, paste0(to, ".var")] <- pred_var
    }

    # output
    return(target)
  }
)

#### log-likelihood ####
setMethod(
  f = "logLik",
  signature = "GP",
  definition = function(object){
    return(object@likelihood)
  }
)

#### Fit ####
# Auxiliary function
.fit_setup <- function(object, contribution, nugget, nugget.t,
                        maxrange, midrange, minrange,
                        azimuth, dip, rake,
                        power, center,
                        contribution.ns,
                        maxrange.ns, midrange.ns, minrange.ns,
                        azimuth.ns, dip.ns, rake.ns){
  # setup
  structures <- sapply(object@model@structures, function(x) x@type)
  Nstruct <- length(structures)
  Ndata <- nrow(object@data)
  data_var <- var(object@data[["value"]])
  data_box <- BoundingBox(object@data)
  data_rbase <- sqrt(sum(data_box[1, ] - data_box[2, ])^2)
  stationary <- class(object@model@structures[[1]]) == "covarianceStructure3D"
  Ncp <- 0 # non stationary control points
  if (!stationary) Ncp <- nrow(object@model@structures[[1]]@contribution@data)

  # parameter vector
  params <- character()
  blocks <- numeric()

  if (contribution.ns & !stationary){
    params <- c(params, rep("contribution.ns.val", Ncp), "contribution.ns.cont",
                "contribution.ns.range", "contribution.ns.mean")
    blocks <- c(blocks, rep(2, Ncp), 1, 1, 1)
  }
  else if (contribution){
    params <- c(params, paste0("contribution.", seq(Nstruct)))
    blocks <- c(blocks, rep(1, length(structures)))
  }
  if (maxrange.ns & !stationary){
    params <- c(params, rep("maxrange.ns.val", Ncp), "maxrange.ns.cont",
                "maxrange.ns.range", "maxrange.ns.mean")
    blocks <- c(blocks, rep(2, Ncp), 1, 1, 1)
  }

  else if (maxrange){
    params <- c(params, paste0("maxrange.", seq(Nstruct)))
    blocks <- c(blocks, rep(1, length(structures)))
  }
  if (midrange.ns & !stationary){
    params <- c(params, rep("midrange.ns.val", Ncp), "midrange.ns.cont",
                "midrange.ns.range", "midrange.ns.mean")
    blocks <- c(blocks, rep(2, Ncp), 1, 1, 1)
  }

  else if (midrange){
    params <- c(params, paste0("midrange.", seq(Nstruct)))
    blocks <- c(blocks, rep(1, length(structures)))
  }
  if (minrange.ns & !stationary){
    params <- c(params, rep("minrange.ns.val", Ncp), "minrange.ns.cont",
                "minrange.ns.range", "minrange.ns.mean")
    blocks <- c(blocks, rep(2, Ncp), 1, 1, 1)
  }
  else if (minrange){
    params <- c(params, paste0("minrange.", seq(Nstruct)))
    blocks <- c(blocks, rep(1, length(structures)))
  }
  if (azimuth.ns & !stationary){
    params <- c(params, rep("azimuth.ns.val", Ncp), "azimuth.ns.cont",
                "azimuth.ns.range", "azimuth.ns.mean")
    blocks <- c(blocks, rep(2, Ncp), 1, 1, 1)
  }
  else if (azimuth){
    params <- c(params, paste0("azimuth.", seq(Nstruct)))
    blocks <- c(blocks, rep(1, length(structures)))
  }
  if (dip.ns & !stationary){
    params <- c(params, rep("dip.ns.val", Ncp), "dip.ns.cont",
                "dip.ns.range", "dip.ns.mean")
    blocks <- c(blocks, rep(2, Ncp), 1, 1, 1)
  }
  else if (dip){
    params <- c(params, paste0("dip.", seq(Nstruct)))
    blocks <- c(blocks, rep(1, length(structures)))
  }
  if (rake.ns & !stationary){
    params <- c(params, rep("rake.ns.val", Ncp), "rake.ns.cont",
                "rake.ns.range", "rake.ns.mean")
    blocks <- c(blocks, rep(2, Ncp), 1, 1, 1)
  }
  else if (rake){
    params <- c(params, paste0("rake.", seq(Nstruct)))
    blocks <- c(blocks, rep(1, length(structures)))
  }

  if (power){
    params <- c(params, paste0("power.", seq(Nstruct)))
    blocks <- c(blocks, rep(1, length(structures)))
  }

  if (center){
    params <- c(params,
                paste0("centerX.", seq(Nstruct)),
                paste0("centerY.", seq(Nstruct)),
                paste0("centerZ.", seq(Nstruct)))
    blocks <- c(blocks, rep(1, length(structures) * 3))
  }

  if (nugget){
    params <- c(params, "nugget")
    blocks <- c(blocks, 1)
  }
  if (nugget.t){
    params <- c(params, "nugget.t")
    blocks <- c(blocks, 1)
  }

  # setting parameters
  opt_min <- opt_max <- xstart <- numeric(length(params))
  # stationary
  if (stationary){
    for (i in seq(Nstruct)) {
      opt_min[params == paste0("contribution.", i)] <- log(data_var / 1000)
      opt_max[params == paste0("contribution.", i)] <- log(data_var * 2)
      xstart[params == paste0("contribution.", i)] <-
        log(object@model@structures[[i]]@params$contribution + 1e-9)

      if (structures[i] != "bias"){
        opt_min[params == paste0("maxrange.", i)] <- data_rbase / 1000
        opt_max[params == paste0("maxrange.", i)] <- data_rbase * 10
        xstart[params == paste0("maxrange.", i)] <-
          object@model@structures[[i]]@params$maxrange

        opt_min[params == paste0("midrange.", i)] <- 0.01
        opt_max[params == paste0("midrange.", i)] <- 1
        xstart[params == paste0("midrange.", i)] <-
          object@model@structures[[i]]@params$midrange /
          object@model@structures[[i]]@params$maxrange

        opt_min[params == paste0("minrange.", i)] <- 0.01
        opt_max[params == paste0("minrange.", i)] <- 1
        xstart[params == paste0("minrange.", i)] <-
          object@model@structures[[i]]@params$minrange /
          object@model@structures[[i]]@params$midrange

        opt_min[params == paste0("azimuth.", i)] <- 0
        opt_max[params == paste0("azimuth.", i)] <- 360
        xstart[params == paste0("azimuth.", i)] <-
          object@model@structures[[i]]@params$azimuth

        opt_min[params == paste0("dip.", i)] <- 0
        opt_max[params == paste0("dip.", i)] <- 90
        xstart[params == paste0("dip.", i)] <-
          object@model@structures[[i]]@params$dip

        opt_min[params == paste0("rake.", i)] <- -90
        opt_max[params == paste0("rake.", i)] <- 90
        xstart[params == paste0("rake.", i)] <-
          object@model@structures[[i]]@params$rake

        if (structures[i] == "cauchy"){
          opt_min[params == paste0("power.", i)] <- 0.1
          opt_max[params == paste0("power.", i)] <- 3
          xstart[params == paste0("power.", i)] <-
            object@model@structures[[i]]@params$power
        }

        if (structures[i] == "nnet"){
          opt_min[params == paste0("power.", i)] <- 0
          opt_max[params == paste0("power.", i)] <- 6
          xstart[params == paste0("power.", i)] <-
            log10(object@model@structures[[i]]@params$power)
        }

        if (structures[i] %in% c("linear", "nnet")){
          opt_min[params == paste0("centerX.", i)] <- data_box[1, 1]
          opt_max[params == paste0("centerX.", i)] <- data_box[2, 1]
          xstart[params == paste0("centerX.", i)] <-
            object@model@structures[[i]]@params$center[1]

          opt_min[params == paste0("centerY.", i)] <- data_box[1, 2]
          opt_max[params == paste0("centerY.", i)] <- data_box[2, 2]
          xstart[params == paste0("centerY.", i)] <-
            object@model@structures[[i]]@params$center[2]

          opt_min[params == paste0("centerZ.", i)] <- data_box[1, 3]
          opt_max[params == paste0("centerZ.", i)] <- data_box[2, 3]
          xstart[params == paste0("centerZ.", i)] <-
            object@model@structures[[i]]@params$center[3]
        }

      }

    }
  }

  opt_min[params == "nugget"] <- log(data_var / 1000)
  opt_max[params == "nugget"] <- log(data_var * 2)
  xstart[params == "nugget"] <- log(object@model@nugget + 1e-9)

  opt_min[params == "nugget.t"] <- data_var / 1000
  opt_max[params == "nugget.t"] <- data_var * 2
  xstart[params == "nugget.t"] <- object@model@nugget.dir

  # non stationary
  opt_min[params == "contribution.ns.val"] <- rep(-5, Ncp)
  opt_max[params == "contribution.ns.val"] <- rep(5, Ncp)
  xstart[params == "contribution.ns.val"] <- rep(0, Ncp)
  opt_min[params == "contribution.ns.cont"] <- log((data_var / 10) ^ 2)
  opt_max[params == "contribution.ns.cont"] <- log((data_var / 4) ^ 2)
  xstart[params == "contribution.ns.cont"] <- log((data_var / 6) ^ 2)
  opt_min[params == "contribution.ns.range"] <- log(data_rbase) - 5
  opt_max[params == "contribution.ns.range"] <- log(data_rbase) + 3
  xstart[params == "contribution.ns.range"] <- log(data_rbase) - 2
  opt_min[params == "contribution.ns.mean"] <- log((data_var / 10) ^ 1)
  opt_max[params == "contribution.ns.mean"] <- log((data_var / 1) ^ 1)
  xstart[params == "contribution.ns.mean"] <- log((data_var / 6) ^ 1)

  opt_min[params == "maxrange.ns.val"] <- rep(-5, Ncp)
  opt_max[params == "maxrange.ns.val"] <- rep(5, Ncp)
  xstart[params == "maxrange.ns.val"] <- rep(0, Ncp)
  opt_min[params == "maxrange.ns.cont"] <- -3
  opt_max[params == "maxrange.ns.cont"] <- 1
  xstart[params == "maxrange.ns.cont"] <- 0
  opt_min[params == "maxrange.ns.range"] <- log(data_rbase) - 5
  opt_max[params == "maxrange.ns.range"] <- log(data_rbase) + 3
  xstart[params == "maxrange.ns.range"] <- log(data_rbase) - 2
  opt_min[params == "maxrange.ns.mean"] <- log(data_rbase) - 5
  opt_max[params == "maxrange.ns.mean"] <- log(data_rbase)
  xstart[params == "maxrange.ns.mean"] <- log(data_rbase) - 3

  opt_min[params == "midrange.ns.val"] <- rep(-5, Ncp)
  opt_max[params == "midrange.ns.val"] <- rep(5, Ncp)
  xstart[params == "midrange.ns.val"] <- rep(0, Ncp)
  opt_min[params == "midrange.ns.cont"] <- 0.1
  opt_max[params == "midrange.ns.cont"] <- 3
  xstart[params == "midrange.ns.cont"] <- 1
  opt_min[params == "midrange.ns.range"] <- log(data_rbase) - 5
  opt_max[params == "midrange.ns.range"] <- log(data_rbase) + 3
  xstart[params == "midrange.ns.range"] <- log(data_rbase) - 2
  opt_min[params == "midrange.ns.mean"] <- -5
  opt_max[params == "midrange.ns.mean"] <- 5
  xstart[params == "midrange.ns.mean"] <- 0

  opt_min[params == "minrange.ns.val"] <- rep(-5, Ncp)
  opt_max[params == "minrange.ns.val"] <- rep(5, Ncp)
  xstart[params == "minrange.ns.val"] <- rep(0, Ncp)
  opt_min[params == "minrange.ns.cont"] <- 0.1
  opt_max[params == "minrange.ns.cont"] <- 3
  xstart[params == "minrange.ns.cont"] <- 1
  opt_min[params == "minrange.ns.range"] <- log(data_rbase) - 5
  opt_max[params == "minrange.ns.range"] <- log(data_rbase) + 3
  xstart[params == "minrange.ns.range"] <- log(data_rbase) - 2
  opt_min[params == "minrange.ns.mean"] <- -5
  opt_max[params == "minrange.ns.mean"] <- 5
  xstart[params == "minrange.ns.mean"] <- 0

  opt_min[params == "azimuth.ns.val"] <- rep(-5, Ncp)
  opt_max[params == "azimuth.ns.val"] <- rep(5, Ncp)
  xstart[params == "azimuth.ns.val"] <- rep(0, Ncp)
  opt_min[params == "azimuth.ns.cont"] <- 0.1
  opt_max[params == "azimuth.ns.cont"] <- 3
  xstart[params == "azimuth.ns.cont"] <- 1
  opt_min[params == "azimuth.ns.range"] <- log(data_rbase) - 5
  opt_max[params == "azimuth.ns.range"] <- log(data_rbase) + 3
  xstart[params == "azimuth.ns.range"] <- log(data_rbase) - 2
  opt_min[params == "azimuth.ns.mean"] <- -5
  opt_max[params == "azimuth.ns.mean"] <- 5
  xstart[params == "azimuth.ns.mean"] <- 0

  opt_min[params == "dip.ns.val"] <- rep(-5, Ncp)
  opt_max[params == "dip.ns.val"] <- rep(5, Ncp)
  xstart[params == "dip.ns.val"] <- rep(0, Ncp)
  opt_min[params == "dip.ns.cont"] <- 0.1
  opt_max[params == "dip.ns.cont"] <- 3
  xstart[params == "dip.ns.cont"] <- 1
  opt_min[params == "dip.ns.range"] <- log(data_rbase) - 5
  opt_max[params == "dip.ns.range"] <- log(data_rbase) + 3
  xstart[params == "dip.ns.range"] <- log(data_rbase) - 2
  opt_min[params == "dip.ns.mean"] <- -5
  opt_max[params == "dip.ns.mean"] <- 5
  xstart[params == "dip.ns.mean"] <- 0

  opt_min[params == "rake.ns.val"] <- rep(-5, Ncp)
  opt_max[params == "rake.ns.val"] <- rep(5, Ncp)
  xstart[params == "rake.ns.val"] <- rep(0, Ncp)
  opt_min[params == "rake.ns.cont"] <- 0.1
  opt_max[params == "rake.ns.cont"] <- 3
  xstart[params == "rake.ns.cont"] <- 1
  opt_min[params == "rake.ns.range"] <- log(data_rbase) - 5
  opt_max[params == "rake.ns.range"] <- log(data_rbase) + 3
  xstart[params == "rake.ns.range"] <- log(data_rbase) - 2
  opt_min[params == "rake.ns.mean"] <- -5
  opt_max[params == "rake.ns.mean"] <- 5
  xstart[params == "rake.ns.mean"] <- 0

  # conforming starting point to limits
  xstart[xstart < opt_min] <- opt_min[xstart < opt_min]
  xstart[xstart > opt_max] <- opt_max[xstart > opt_max]

  # output
  return(list(params = params, opt_min = opt_min, opt_max = opt_max,
              xstart = xstart, blocks = blocks,
              stationary = stationary))
}

.build_cov <- function(x, optdata, model){
  m <- model@structures

  if (optdata$stationary){
    # covariance structures
    for(i in seq_along(m)){
      if (paste0("contribution.", i) %in% optdata$params){
        pos <- which(optdata$params == paste0("contribution.", i))
        m[[i]]@params$contribution <- exp(x[pos])
      }
      if (m[[i]]@type != "bias"){
        if (paste0("maxrange.", i) %in% optdata$params){
          pos <- which(optdata$params == paste0("maxrange.", i))
          m[[i]]@params$maxrange <- x[pos]
        }
        if (paste0("midrange.", i) %in% optdata$params){
          pos <- which(optdata$params == paste0("midrange.", i))
          m[[i]]@params$midrange <- x[pos] * m[[i]]@params$maxrange
        }
        else
          m[[i]]@params$midrange <- m[[i]]@params$maxrange
        if (paste0("minrange.", i) %in% optdata$params){
          pos <- which(optdata$params == paste0("minrange.", i))
          m[[i]]@params$minrange <- x[pos] * m[[i]]@params$midrange
        }
        else
          m[[i]]@params$minrange <- m[[i]]@params$midrange
        if (paste0("azimuth.", i) %in% optdata$params){
          pos <- which(optdata$params == paste0("azimuth.", i))
          m[[i]]@params$azimuth <- x[pos]
        }
        if (paste0("dip.", i) %in% optdata$params){
          pos <- which(optdata$params == paste0("dip.", i))
          m[[i]]@params$dip <- x[pos]
        }
        if (paste0("rake.", i) %in% optdata$params){
          pos <- which(optdata$params == paste0("rake.", i))
          m[[i]]@params$rake <- x[pos]
        }

        if (m[[i]]@type == "cauchy"){
          if (paste0("power.", i) %in% optdata$params){
            pos <- which(optdata$params == paste0("power.", i))
            m[[i]]@params$power <- x[pos]
          }
        }

        if (m[[i]]@type == "nnet"){
          if (paste0("power.", i) %in% optdata$params){
            pos <- which(optdata$params == paste0("power.", i))
            m[[i]]@params$power <- 10 ^ x[pos]
          }
        }

        if (m[[i]]@type %in% c("linear", "nnet")){
          if (paste0("centerX.", i) %in% optdata$params){
            pos <- which(optdata$params %in% c(paste0("centerX.", i),
                                               paste0("centerY.", i),
                                               paste0("centerZ.", i)))
            m[[i]]@params$center <- x[pos]
          }
        }
      }
    }
  }
  else{
    # type <- m[[1]]@type
    type <- "gaussian"
    reg <- 1e-3

    # contribution
    tmpdata <- m[[1]]@contribution@data
    if ("contribution.ns.mean" %in% optdata$params){
      avg <- x[which(optdata$params == "contribution.ns.mean")]
      cont <- x[which(optdata$params == "contribution.ns.cont")]
      rng <- x[which(optdata$params == "contribution.ns.range")]
      z <- x[which(optdata$params == "contribution.ns.val")]
      tmpgp <- SPGP(data = points3DDataFrame(),
                  model = covarianceModel3D(
                    0, covarianceStructure3D(type, exp(cont), exp(rng))),
                  mean = avg, reg.v = reg,
                  pseudo_inputs = m[[1]]@contribution@pseudo_inputs)
      tmpdata <- Simulate(tmpgp, tmpdata, to = "value", Nsim = 1,
                          discount.noise = T, verbose = F, randnum = z)
      m[[1]]@contribution <- SPGP(data = tmpdata, value = "value.sim_1",
                                  model = covarianceModel3D(
                                    0, covarianceStructure3D(type, exp(cont), exp(rng))),
                                  mean = avg, reg.v = reg,
                                  pseudo_inputs = m[[1]]@contribution@pseudo_inputs)
    }
    else if ("contribution.1" %in% optdata$params){
      tmpdata[, "value"] <- x[which(optdata$params == "contribution.1")]
      m[[1]]@contribution <- SPGP(data = tmpdata, value = "value",
                                  model = covarianceModel3D(
                                    0, covarianceStructure3D(type, 0, 1)),
                                  reg.v = reg)
    }

    # maxrange
    tmpdata <- m[[1]]@maxrange@data
    if ("maxrange.ns.mean" %in% optdata$params){
      avg <- x[which(optdata$params == "maxrange.ns.mean")]
      cont <- x[which(optdata$params == "maxrange.ns.cont")]
      rng <- x[which(optdata$params == "maxrange.ns.range")]
      z <- x[which(optdata$params == "maxrange.ns.val")]
      tmpgp <- SPGP(data = points3DDataFrame(),
                    model = covarianceModel3D(
                      0, covarianceStructure3D(type, exp(cont), exp(rng))),
                    mean = avg, reg.v = reg,
                    pseudo_inputs = m[[1]]@maxrange@pseudo_inputs)
      tmpdata <- Simulate(tmpgp, tmpdata, to = "value", Nsim = 1,
                          discount.noise = T, verbose = F, randnum = z)
      # if (any(is.nan(tmpdata[["value"]]))) stop("Não deu")
      m[[1]]@maxrange <- SPGP(data = tmpdata, value = "value.sim_1",
                                  model = covarianceModel3D(
                                    0, covarianceStructure3D(type, exp(cont), exp(rng))),
                                  mean = avg, reg.v = reg,
                              pseudo_inputs = m[[1]]@contribution@pseudo_inputs)
    }
    else if ("maxrange.1" %in% optdata$params){
      tmpdata[, "value"] <- x[which(optdata$params == "maxrange.1")]
      m[[1]]@maxrange <- SPGP(data = tmpdata, value = "value",
                                  model = covarianceModel3D(
                                    0, covarianceStructure3D(type, 0, 1)),
                                  reg.v = reg)
    }

    # midrange
    tmpdata <- m[[1]]@midrange@data
    if ("midrange.ns.mean" %in% optdata$params){
      avg <- x[which(optdata$params == "midrange.ns.mean")]
      cont <- x[which(optdata$params == "midrange.ns.cont")]
      rng <- x[which(optdata$params == "midrange.ns.range")]
      z <- x[which(optdata$params == "midrange.ns.val")]
      tmpgp <- SPGP(data = points3DDataFrame(),
                    model = covarianceModel3D(
                      0, covarianceStructure3D(type, cont, exp(rng))),
                    mean = avg, reg.v = reg,
                    pseudo_inputs = m[[1]]@midrange@pseudo_inputs)
      tmpdata <- Simulate(tmpgp, tmpdata, to = "value", Nsim = 1,
                          discount.noise = T, verbose = F, randnum = z)
      # if (any(is.nan(tmpdata[["value.sim_1"]]))) browser()
      m[[1]]@midrange <- SPGP(data = tmpdata, value = "value.sim_1",
                              model = covarianceModel3D(
                                0, covarianceStructure3D(type, cont, exp(rng))),
                              mean = avg, reg.v = reg,
                              pseudo_inputs = m[[1]]@contribution@pseudo_inputs)
    }
    else if ("midrange.1" %in% optdata$params){
      tmpdata[, "value"] <- x[which(optdata$params == "midrange.1")]
      m[[1]]@midrange <- SPGP(data = tmpdata, value = "value",
                              model = covarianceModel3D(
                                0, covarianceStructure3D(type, 0, 1)),
                              reg.v = reg)
    }

    # minrange
    tmpdata <- m[[1]]@minrange@data
    if ("minrange.ns.mean" %in% optdata$params){
      avg <- x[which(optdata$params == "minrange.ns.mean")]
      cont <- x[which(optdata$params == "minrange.ns.cont")]
      rng <- x[which(optdata$params == "minrange.ns.range")]
      z <- x[which(optdata$params == "minrange.ns.val")]
      tmpgp <- SPGP(data = points3DDataFrame(),
                    model = covarianceModel3D(
                      0, covarianceStructure3D(type, cont, exp(rng))),
                    mean = avg, reg.v = reg,
                    pseudo_inputs = m[[1]]@minrange@pseudo_inputs)
      tmpdata <- Simulate(tmpgp, tmpdata, to = "value", Nsim = 1,
                          discount.noise = T, verbose = F, randnum = z)
      m[[1]]@minrange <- SPGP(data = tmpdata, value = "value.sim_1",
                              model = covarianceModel3D(
                                0, covarianceStructure3D(type, cont, exp(rng))),
                              mean = avg, reg.v = reg,
                              pseudo_inputs = m[[1]]@contribution@pseudo_inputs)
    }
    else if ("minrange.1" %in% optdata$params){
      tmpdata[, "value"] <- x[which(optdata$params == "minrange.1")]
      m[[1]]@minrange <- SPGP(data = tmpdata, value = "value",
                              model = covarianceModel3D(
                                0, covarianceStructure3D(type, 0, 1)),
                              reg.v = reg)
    }

    # azimuth
    tmpdata <- m[[1]]@azimuth@data
    if ("azimuth.ns.mean" %in% optdata$params){
      avg <- x[which(optdata$params == "azimuth.ns.mean")]
      cont <- x[which(optdata$params == "azimuth.ns.cont")]
      rng <- x[which(optdata$params == "azimuth.ns.range")]
      z <- x[which(optdata$params == "azimuth.ns.val")]
      tmpgp <- SPGP(data = points3DDataFrame(),
                    model = covarianceModel3D(
                      0, covarianceStructure3D(type, cont, exp(rng))),
                    mean = avg, reg.v = reg,
                    pseudo_inputs = m[[1]]@azimuth@pseudo_inputs)
      tmpdata <- Simulate(tmpgp, tmpdata, to = "value", Nsim = 1,
                          discount.noise = T, verbose = F, randnum = z)
      # if (any(is.nan(tmpdata[["value.sim_1"]]))) browser()
      m[[1]]@azimuth <- SPGP(data = tmpdata, value = "value.sim_1",
                              model = covarianceModel3D(
                                0, covarianceStructure3D(type, cont, exp(rng))),
                              mean = avg, reg.v = reg,
                             pseudo_inputs = m[[1]]@contribution@pseudo_inputs)
    }
    else if ("azimuth.1" %in% optdata$params){
      tmpdata[, "value"] <- x[which(optdata$params == "azimuth.1")]
      m[[1]]@azimuth <- SPGP(data = tmpdata, value = "value",
                              model = covarianceModel3D(
                                0, covarianceStructure3D(type, 0, 1)),
                              reg.v = reg)
    }

    # dip
    tmpdata <- m[[1]]@dip@data
    if ("dip.ns.mean" %in% optdata$params){
      avg <- x[which(optdata$params == "dip.ns.mean")]
      cont <- x[which(optdata$params == "dip.ns.cont")]
      rng <- x[which(optdata$params == "dip.ns.range")]
      z <- x[which(optdata$params == "dip.ns.val")]
      tmpgp <- SPGP(data = points3DDataFrame(),
                    model = covarianceModel3D(
                      0, covarianceStructure3D(type, cont, exp(rng))),
                    mean = avg, reg.v = reg,
                    pseudo_inputs = m[[1]]@dip@pseudo_inputs)
      tmpdata <- Simulate(tmpgp, tmpdata, to = "value", Nsim = 1,
                          discount.noise = T, verbose = F, randnum = z)
      m[[1]]@dip <- SPGP(data = tmpdata, value = "value.sim_1",
                              model = covarianceModel3D(
                                0, covarianceStructure3D(type, cont, exp(rng))),
                              mean = avg, reg.v = reg,
                         pseudo_inputs = m[[1]]@contribution@pseudo_inputs)
    }
    else if ("dip.1" %in% optdata$params){
      tmpdata[, "value"] <- x[which(optdata$params == "dip.1")]
      m[[1]]@dip <- SPGP(data = tmpdata, value = "value",
                              model = covarianceModel3D(
                                0, covarianceStructure3D(type, 0, 1)),
                              reg.v = reg)
    }

    # rake
    tmpdata <- m[[1]]@rake@data
    if ("rake.ns.mean" %in% optdata$params){
      avg <- x[which(optdata$params == "rake.ns.mean")]
      cont <- x[which(optdata$params == "rake.ns.cont")]
      rng <- x[which(optdata$params == "rake.ns.range")]
      z <- x[which(optdata$params == "rake.ns.val")]
      tmpgp <- SPGP(data = points3DDataFrame(),
                    model = covarianceModel3D(
                      0, covarianceStructure3D(type, cont, exp(rng))),
                    mean = avg, reg.v = reg,
                    pseudo_inputs = m[[1]]@rake@pseudo_inputs)
      tmpdata <- Simulate(tmpgp, tmpdata, to = "value", Nsim = 1,
                          discount.noise = T, verbose = F, randnum = z)
      m[[1]]@rake <- SPGP(data = tmpdata, value = "value.sim_1",
                              model = covarianceModel3D(
                                0, covarianceStructure3D(type, cont, exp(rng))),
                              mean = avg, reg.v = reg,
                          pseudo_inputs = m[[1]]@contribution@pseudo_inputs)
    }
    else if ("rake.1" %in% optdata$params){
      tmpdata[, "value"] <- x[which(optdata$params == "rake.1")]
      m[[1]]@rake <- SPGP(data = tmpdata, value = "value",
                              model = covarianceModel3D(
                                0, covarianceStructure3D(type, 0, 1)),
                              reg.v = reg)
    }
  }


  # nugget
  if ("nugget" %in% optdata$params){
    pos <- which(optdata$params == "nugget")
    tmpnug <- exp(x[pos])
  }
  else
    tmpnug <- model@nugget


  # nugget for tangents
  if ("nugget.t" %in% optdata$params){
    pos <- which(optdata$params == "nugget.t")
    tmpnug.t <- x[pos]
  }
  else
    tmpnug.t <- model@nugget.dir

  # output
  return(covarianceModel3D(tmpnug, m, tmpnug.t))
}


#' @rdname Fit
setMethod(
  f = "Fit",
  signature = "GP",
  definition = function(object, contribution = T, nugget = T, nugget.t = F,
                        maxrange = T, midrange = F, minrange = F,
                        azimuth = F, dip = F, rake = F,
                        power = F, center = F,
                        contribution.ns = F,
                        maxrange.ns = F, midrange.ns = F, minrange.ns = F,
                        azimuth.ns = F, dip.ns = F, rake.ns = F, ...){

    # setup
    int <- object@data[["interpolate"]]
    optdata <- .fit_setup(object, contribution, nugget, nugget.t,
                          maxrange, midrange, minrange,
                          azimuth, dip, rake,
                          power, center,
                          contribution.ns,
                          maxrange.ns, midrange.ns, minrange.ns,
                          azimuth.ns, dip.ns, rake.ns)

    # fitness function
    makeGP <- function(x, finished = F){
      model <- .build_cov(x, optdata, object@model)

      tmpgp <- GP(
        data = object@data,
        model = model,
        value = "value",
        mean = object@mean,
        trend = object@trend,
        tangents = object@tangents,
        force.interp = int
      )
      # output
      if(finished)
        return(tmpgp)
      else
        return(tmpgp@likelihood) # + logLik(model)
    }

    # optimization
    opt <- GeneticTrainingReal(
      fitness = function(x) makeGP(x, F),
      minval = optdata$opt_min,
      maxval = optdata$opt_max,
      start = optdata$xstart,
      blocks = optdata$blocks,
      ...
    )

    # update
    sol <- opt$bestsol
    return(makeGP(sol, T))
  }
)

#### Simulate ####
#' @rdname Simulate
setMethod(
  f = "Simulate",
  signature = "GP",
  definition = function(object, target, to = "value", Nsim,
                        discount.noise = F){

    # setup
    w_var <- object@pre_comp$w_var
    w_value <- object@pre_comp$w_value

    # number formatting
    Ndigits <- length(unlist(strsplit(as.character(Nsim), NULL)))
    form <- paste0("%0", Ndigits, ".0f")

    # trivial solution
    # if(object@model@total.var == 0){
    #   sim <- matrix(object@mean, nrow(target), Nsim)
    #   colnames(sim) <- paste0(to, ".sim_", sprintf(form, seq(Nsim)))
    #   target[, colnames(sim)] <- data.frame(as.matrix(sim))
    #   return(target)
    # }

    # covariances
    Ntang <- nrow(object@tangents)
    Ktarget <- CovarianceMatrix(target, object@data,
                                object@model)
    if(Ntang > 0){
      K1 <- CovarianceMatrixD1(target,
                               object@tangents,
                               object@model)
      Ktarget <- cbind(Ktarget, K1)
    }
    # pre-computation
    LinvK <- solve(w_var, t(Ktarget))

    # prediction
    # residuals
    pred <- apply(Ktarget, 1, function(rw)sum(rw * w_value)) + object@mean
    # trend
    if(length(object@trend) > 0){
      TRtarget <- TrendMatrix(target, object@trend)
      TRdata <- TrendMatrix(object@data, object@trend)
      w_tr <- object@pre_comp$w_trend
      beta <- object@beta
      R <- t(TRtarget) - w_tr %*% LinvK
      pred <- pred + t(R) %*% beta
    }

    # covariance matrix
    cov_reduction <- crossprod(LinvK)
    prior_cov <- CovarianceMatrix(target, target, object@model) #+
      #diag(1e-6, nrow(target), nrow(target)) # regularization
    if(!discount.noise)
      prior_cov <- prior_cov + diag(object@model@nugget, nrow(target), nrow(target))
    pred_cov <- prior_cov - cov_reduction
    if(length(object@trend) > 0){
      A <- object@pre_comp$A
      tr_cov <- t(R) %*% solve(A, R)
      tr_cov <- 0.5 * tr_cov + 0.5 * t(tr_cov) # enforcing symmetry
      pred_cov <- pred_cov + tr_cov
    }
    # pred_cov <- 0.5 * pred_cov + 0.5 * t(pred_cov)
    Lpred <- .safeChol(pred_cov)

    # simulation
    rand <- matrix(rnorm(nrow(target) * Nsim), nrow(target), Nsim)
    sim <- Lpred %*% rand + matrix(pred, nrow(target), Nsim)
    colnames(sim) <- paste0(to, ".sim_", sprintf(form, seq(Nsim)))
    target[, colnames(sim)] <- data.frame(as.matrix(sim))

    return(target)
  }
)

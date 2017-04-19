#' @include generics.R
#' @include utils.R
#' @include spatial3DDataFrame.R
#' @include points3DDataFrame.R
NULL

#### Gaussian Process class ####
#' Gaussian Process
#'
#' Implementation of the Gaussian Process model for 3D spatial interpolation.
#'
#' @slot data A \code{spatial3DDataFrame} object containing the necessary data.
#' @slot tangents A \code{points3DDataFrame} object containing structural
#' geology data. Most likely generated with the \code{GetLineDirections()}
#' method.
#' @slot model The covariance model. A \code{list} containing one or more
#' \code{covarianceStructure3D} objects.
#' @slot nugget The model's nugget effect or noise variance.
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
#' @seealso \code{\link{GP-init}}, \code{\link{GP_geomod-class}}
#'
#' @export GP
GP <- setClass(
  "GP",
  slots = c(data = "spatial3DDataFrame",
            tangents = "points3DDataFrame",
            model = "list",
            nugget = "numeric",
            mean = "numeric",
            trend = "character",
            beta = "matrix",
            likelihood = "numeric",
            pre_comp = "list"),
  validity = function(object) {
    if(!all(rapply(object@model,class) == "covarianceStructure3D"))
      stop("Invalid covariance object")
  }
)

#### initialization ####
#' Gaussian Process
#'
#' Implementation of the Gaussian Process model for 3D spatial interpolation.
#'
#' @param data A \code{spatial3DDataFrame} object containing the data one
#' wishes to model.
#' @param model The covariance model. A \code{covarianceStructure3D} or a
#' \code{list} containing one or more such objects.
#' @param value The column name of the variable to be modeled. It is assumed
#' the column does not contain missing values.
#' @param nugget The model's nugget effect or noise variance.
#' @param mean The global mean. Irrelevant if a trend is provided.
#' @param trend The model's trend component. A formula in character format.
#' @param weights The importance of each data point in the model (a vector with
#' values between 0 and 1)
#' @param force.interp Indices of points that must be interpolated exactly.
#' @param reg.v Regularization to improve stability. A single value or a vector
#' with length matching the number of data points.
#' @param tangents A \code{points3DDataFrame} object containing structural
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
#' @seealso \code{\link{GP-class}}
setMethod(
  f = "initialize",
  signature = "GP",
  definition = function(.Object, data, model, value, nugget,
                        mean = NULL, trend = NULL, weights = NULL,
                        force.interp = numeric(), reg.v = 1e-9,
                        tangents = NULL, reg.t = 1e-12){

    # covariance model
    if(length(model) == 1 & class(model) != "list")
      .Object@model <- list(model)
    else
      .Object@model <- model

    # value
    if(length(value) == 1 & class(value) == "character")
      data["value"] <- data[value]
    else
      data["value"] <- value
    yval <- data[["value"]]

    # nugget
    .Object@nugget <- nugget

    # weights
    if(is.null(weights)) weights <- rep(1, nrow(data))
    data["weights"] <- weights

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
    .Object@data <- data[c("value", "weights", "interpolate")]
    .Object@tangents <- as(tangents, "points3DDataFrame")

    # trend
    if(is.null(trend) | length(trend) == 0){
      TR <- matrix(0,nrow(data),0)
    }else{
      TR <- TrendMatrix(data, trend)
      if(!is.null(tangents) && nrow(tangents) > 0){
        TR <- rbind(
          TR,
          TrendMatrixD1(tangents, trend)
        )
      }
      mean <- 0
    }
    Ntrend <- dim(TR)[2]

    # covariances
    Ntang <- 0
    nug2 <-  nugget / (weights + reg.v)
    nug2[force.interp] <- reg.v
    K <- CovarianceMatrix(data, data, model, T) +
      diag(nug2, length(nug2), length(nug2))
    if(!is.null(tangents) && nrow(tangents) > 0){
      K1 <- CovarianceMatrixD1(data, tangents, model, T)
      K2 <- CovarianceMatrixD2(tangents, model, T)
      Ntang <- nrow(K2)
      # regularization
      K2 <- K2 + diag(reg.t, Ntang, Ntang)
      # final matrix
      K <- rbind(
        cbind(K, K1),
        cbind(t(K1), K2)
      )
    }
    # enforcing symmetry
    K <- 0.5 * K + 0.5 * t(K)

    # pre-computations
    .Object@pre_comp <- list()
    yval <- c(yval - mean, rep(0, Ntang))
    .Object@mean <- mean

    # Rcpp
    # pc <- GP_precomp(K, yval)
    # L <- pc$L
    # LiY <- pc$LiY
    # w_value <- pc$a

    # no Rcpp
    L <- t(chol(Matrix(K)))  # makes L lower triangular
    LiY <- solve(L, Matrix(yval, length(yval), 1))
    w_value <- solve(t(L), LiY)

    .Object@pre_comp$w_value <- as.numeric(w_value)
    .Object@pre_comp$w_var <- L
    if(Ntrend > 0){
      HLi <- solve(L, TR)
      A <- t(HLi) %*% HLi
      b1 <- t(TR) %*% w_value
      beta <- solve(A, b1)
      .Object@beta <- as.matrix(beta)
      .Object@trend <- trend
      .Object@pre_comp$w_trend <- HLi
    }

    # likelihood
    dt <- 2*sum(diag(L)^2)
    .Object@likelihood <- - 0.5 * dt -
      0.5 * sum(yval * w_value) -
      0.5 * length(yval) * log(2 * pi)
    if(Ntrend > 0){
      LiYH <- t(LiY) %*% HLi
      tmp <- LiYH %*% solve(A, t(LiYH))
      .Object@likelihood <- .Object@likelihood +
        0.5 * as.numeric(tmp) +
        Ntrend * log(2*pi) -
        0.5 * determinant(A)$modulus
    }

    # end
    validObject(.Object)
    return(.Object)
  }
)

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

#### Predict ####
#' @rdname Predict
setMethod(
  f = "Predict",
  signature = "GP",
  definition = function(object, target, to = "value", output.var = F){

    # pre processing
    w_var <- object@pre_comp$w_var
    w_value <- object@pre_comp$w_value

    # trend
    if(length(object@trend) > 0){
      TRtarget <- TrendMatrix(target, object@trend)
      TRdata <- TrendMatrix(object@data, object@trend)
      w_tr <- object@pre_comp$w_trend
      beta <- object@beta
    }

    # slicing target to save memory
    maxgrid <- 1000 # optimize this
    Ngrid <- nrow(target)
    Nslice <- ceiling(Ngrid / maxgrid)
    t2 <- Pointify(target)

    for(i in seq(Nslice)){

      # slice ID
      slid <- seq((i - 1) * maxgrid + 1, min(Ngrid, i * maxgrid))
      ttemp <- t2[slid, ]

      # covariances
      Ntang <- nrow(object@tangents)
      Ktarget <- CovarianceMatrix(ttemp, object@data,
                                   object@model, T)
      if(Ntang > 0){
        K1 <- CovarianceMatrixD1(ttemp,
                                   object@tangents,
                                   object@model, T)
        Ktarget <- cbind(Ktarget, K1)
      }

      # prediction
      # residuals
      pred <- apply(Ktarget, 1, function(rw){
        sum(rw * w_value)
      }) + object@mean
      # trend
      if(length(object@trend) > 0){
        LinvK <- solve(w_var, t(Ktarget))
        R <- t(TRtarget[slid,]) - t(w_tr) %*% LinvK
        pred <- pred + t(R) %*% beta
      }
      target[slid, to] <- pred

      # variance
      if(output.var){
        tot_var <- sum(sapply(object@model,
                              function(m) m@contribution))
        if(length(object@trend) > 0){
          pred_var <- colSums(LinvK^2)
          pred_var[pred_var > tot_var] <- tot_var
          pred_var <- tot_var - pred_var + object@nugget
          tr_var <- colSums(
            R * (solve(t(w_tr) %*% w_tr, R))
          )
          pred_var <- pred_var + tr_var
        }
        else{
          pred_var <- colSums(solve(w_var, t(Ktarget))^2)
          pred_var[pred_var > tot_var] <- tot_var
          pred_var <- tot_var - pred_var + object@nugget
        }
        target[slid, paste0(to, ".var")] <- pred_var
      }

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
#' @rdname Fit
setMethod(
  f = "Fit",
  signature = "GP",
  definition = function(object,
                        midrange = F, minrange = F,
                        azimuth = F, dip = F, rake = F,
                        nugget = F, power = F){

    # setup
    structures <- sapply(object@model, function(x) x@type)
    Nstruct <- length(structures)
    Ndata <- nrow(object@data)
    data_var <- var(object@data[["value"]])
    data_box <- BoundingBox(object@data)
    data_rbase <- sqrt(sum(data_box[1, ] - data_box[2, ])^2)
    int <- object@data[["interpolate"]]
    w <- object@data[["weights"]]

    # optimization limits and starting point
    opt_min <- opt_max <- numeric(Nstruct * 8 + 1)
    xstart <- matrix(0, 1, Nstruct * 8 + 1)
    for (i in 1:Nstruct){
      # contribution
      opt_min[(i - 1) * 8 + 1] <- data_var / 1000
      opt_max[(i - 1) * 8 + 1] <- data_var * 2
      xstart[(i - 1) * 8 + 1] <- object@model[[i]]@contribution

      # maxrange
      opt_min[(i - 1) * 8 + 2] <- data_rbase / 1000
      opt_max[(i - 1) * 8 + 2] <- data_rbase
      xstart[(i - 1) * 8 + 2] <- object@model[[i]]@maxrange

      # midrange (multiple of maxrange)
      if (midrange)
        opt_min[(i - 1) * 8 + 3] <- 0.01
      else
        opt_min[(i - 1) * 8 + 3] <- 1
      opt_max[(i - 1) * 8 + 3] <- 1
      xstart[(i - 1) * 8 + 3] <- object@model[[i]]@midrange /
        object@model[[i]]@maxrange

      # minrange(multiple of midrange)
      if (minrange)
        opt_min[(i - 1) * 8 + 4] <- 0.01
      else
        opt_min[(i - 1) * 8 + 4] <- 1
      opt_max[(i - 1) * 8 + 4] <- 1
      xstart[(i - 1) * 8 + 4] <- object@model[[i]]@minrange /
        object@model[[i]]@midrange

      # azimuth
      opt_min[(i - 1) * 8 + 5] <- 0
      if (azimuth)
        opt_max[(i - 1) * 8 + 5] <- 360
      else
        opt_max[(i - 1) * 8 + 5] <- 0
      xstart[(i - 1) * 8 + 5] <- object@model[[i]]@azimuth

      # dip
      opt_min[(i - 1) * 8 + 6] <- 0
      if (dip)
        opt_max[(i - 1) * 8 + 6] <- 90
      else
        opt_max[(i - 1) * 8 + 6] <- 0
      xstart[(i - 1) * 8 + 6] <- object@model[[i]]@dip

      # rake
      opt_min[(i - 1) * 8 + 7] <- 0
      if (rake)
        opt_max[(i - 1) * 8 + 7] <- 90
      else
        opt_max[(i - 1) * 8 + 7] <- 0
      xstart[(i - 1) * 8 + 7] <- object@model[[i]]@rake

      # power
      if (power){
        opt_min[(i - 1) * 8 + 8] <- 0.1
        opt_max[(i - 1) * 8 + 8] <- 3
      }
      else{
        opt_min[(i - 1) * 8 + 8] <- 1
        opt_max[(i - 1) * 8 + 8] <- 1
      }
      xstart[(i - 1) * 8 + 8] <- object@model[[i]]@power
    }

    # nugget
    if (nugget){
      opt_min[Nstruct * 8 + 1] <- data_var / 1000
      opt_max[Nstruct * 8 + 1] <- data_var * 2
    }
    else{
      opt_min[Nstruct * 8 + 1] <- 0 # not used
      opt_max[Nstruct * 8 + 1] <- 0 # not used
    }
    xstart[Nstruct * 8 + 1] <- object@nugget

    # conforming starting point to limits
    xstart[xstart < opt_min] <- opt_min[xstart < opt_min]
    xstart[xstart > opt_max] <- opt_max[xstart > opt_max]

    # fitness function
    makeGP <- function(x, finished = F){
      # covariance model
      m <- vector("list", Nstruct)
      for(i in 1:Nstruct){
        m[[i]] <- covarianceStructure3D(
          type = structures[i],
          contribution = x[(i - 1) * 8 + 1],
          maxrange = x[(i - 1) * 8 + 2],
          midrange = x[(i - 1) * 8 + 2] * x[(i - 1) * 8 + 3],
          minrange = x[(i - 1) * 8 + 2] * x[(i - 1) * 8 + 3] *
            x[(i - 1) * 8 + 4],
          azimuth = x[(i - 1) * 8 + 5],
          dip = x[(i - 1) * 8 + 6],
          rake = x[(i - 1) * 8 + 7],
          power = x[(i - 1) * 8 + 8]
        )
      }
      # nugget
      if(nugget)
              tmpnug <- x[Nstruct * 8 + 1]
      else
              tmpnug <- object@nugget
      # GP
      tmpgp <- GP(
        data = object@data,
        model = m,
        value = "value",
        nugget = tmpnug,
        mean = object@mean,
        trend = object@trend,
        tangents = object@tangents,
        weights = w,
        force.interp = int
      )
      # output
      if(finished)
        return(tmpgp)
      else
        return(tmpgp@likelihood)
    }

    # optimization
    opt <- ga(
      type = "real-valued",
      fitness = function(x) makeGP(x, F),
      min = opt_min,
      max = opt_max,
      pmutation = 0.5,
      popSize = 20,
      run = 20,
      monitor = F,
      suggestions = xstart
    )

    # update
    sol <- opt@solution
    return(makeGP(sol, T))
  }
)

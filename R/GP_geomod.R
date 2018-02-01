#' @include generics.R
#' @include utils.R
#' @include spatial3DDataFrame.R
#' @include points3DDataFrame.R
NULL

#' Gaussian Process for Implicit Geological Modeling
#'
#' Implicit geological using Gaussian Processes and compositional relations
#' among the geological classes.
#'
#' @slot GPs A \code{list} containing one \code{GP} object per class.
#' @slot params A \code{list} containing modeling parameters such as
#' regularization factors.
#'
#' @references Gonçalves IG, Kumaira S, Guadagnin F.
#' A machine learning approach to the potential-field method for implicit
#' modeling of geological structures. Comput Geosci 2017;103:173–82.
#' doi:10.1016/j.cageo.2017.03.015.
#'
#' @seealso \code{\link{GP_geomod-init}}, \code{\link{GP-class}}
#'
#' @name GP_geomod-class
#' @export GP_geomod
GP_geomod <- setClass(
  "GP_geomod",
  slots = c(GPs = "list",
            params = "list")
)

#### initialization ####
#' Gaussian Process for Implicit Geological Modeling
#'
#' Implicit geological using Gaussian Processes and compositional relations
#' among the geological classes.
#'
#' @param data A 3D spatial object.
#' @param value1 A column name indicating the variable that contains the
#' geological classes.
#' @param value2 Another column name, with possibly different class labels.
#' @param model The covariance model. A \code{covarianceStructure3D} or a
#' \code{list} containing one or more such objects.
#' @param nugget The model's nugget effect or noise variance.
#' @param tangents A \code{points3DDataFrame} object containing structural
#' geology data. Most likely generated with the \code{GetLineDirections()}
#' method.
#' @param enforce.contacts Force the model to interpolate geological contacts?
#' @param reg.v Regularization to improve stability. A single value or a vector
#' with length matching the number of data points.
#' @param reg.t Regularization for structural data. A single value or a vector
#' with length matching the number of structural data.
#'
#' @details The role of this class is to manage multiple \code{GP} objects, one
#' per geological class in the data, that model an indicator variable or
#' "potential" for each class. A point in 3D space is assigned to the class
#' with the highest modeled potential. The class labels are defined as
#' \code{unique(c(data[[value1]], data[[value2]]))}.
#'
#' The points at which \code{data[[value1]] != data[[value2]]} are considered
#' to lie exactly at a boundary between two geological classes. For these points
#' the \code{force.interp} option of the \code{GP} class is turned on in order
#' to make the geological boundaries pass through them.
#'
#' The prediction provides an indicator value for each class at each location,
#' as well as the label of the most likely class. In locations away from the
#' data, the most likely class can be an extra, \code{"Unknown"} class. This
#' class is there to to ensure that the probabilities sum to 1 according to
#' compositional data principles, and also as a measure of model uncertainty.
#'
#' The geological boundaries can be drawn by contouring each class's indicator
#' at the value 0.
#'
#' @seealso \code{\link{GP}}, \code{\link{GP-init}}, \code{\link{Predict}},
#' \code{\link{Make3DArray}}
#'
#' @name GP_geomod-init
#'
#' @references Gonçalves IG, Kumaira S, Guadagnin F.
#' A machine learning approach to the potential-field method for implicit
#' modeling of geological structures. Comput Geosci 2017;103:173–82.
#' doi:10.1016/j.cageo.2017.03.015.
GP_geomod <- function(data, value1, value2 = value1,
                      model, tangents = NULL,
                      enforce.contacts = T,
                      reg.v = 1e-9, reg.t = 1e-9){
  # setup
  Ndata <- nrow(data)
  xdata <- GetData(data)
  categories <- c(xdata[, value1], xdata[, value2])
  categories[!is.na(categories)]
  categories <- sort(unique(categories))
  ncat <- length(categories)
  GPs <- vector("list", ncat)

  # model building
  for (i in seq_along(categories)){

    # indicators, or potential components
    ind <- matrix(- 1 / ncat, nrow(data), 2) # negative
    ind[xdata[, value1] == categories[i], 1] <- 1 # positive 1
    ind[xdata[, value2] == categories[i], 2] <- 1 # positive 2
    ind <- rowMeans(ind) # contacts get (1 - 1 / ncat) / 2

    # contact indices
    cont <- ind == (1 - 1 / ncat) / 2
    if (!enforce.contacts) cont <- numeric()

    # GPs
    GPs[[i]] <- GP(data, model,
                   ind,
                   mean = - 1 / ncat,
                   tangents = tangents,
                   reg.t = reg.t,
                   reg.v = reg.v,
                   force.interp = cont)
  }
  names(GPs) <- make.names(categories)

  # output
  new("GP_geomod", GPs = GPs,
      params = list(reg.v = reg.v, reg.t = reg.t))
}


#### show ####
setMethod(
  f = "show",
  signature = "GP_geomod",
  definition = function(object){
    # display
    cat("Object of class ", class(object), "\n", sep = "")
    cat("Data points:", nrow(object@GPs[[1]]@data), "\n")
    Ntang <- ifelse(is.null(object@GPs[[1]]@tangents), 0,
                    nrow(object@GPs[[1]]@tangents))
    cat("Tangent points:", Ntang, "\n")
    cat("Number of classes:", length(object@GPs), "\n")
    cat("Log-likelihood:", logLik(object), "\n")
    lab <- names(object@GPs)
    cat("Class labels:\n")
    for (i in seq_along(lab))
      cat("  ", lab[i], "\n")
  }
)

#### log-likelihood ####
setMethod(
  f = "logLik",
  signature = "GP_geomod",
  definition = function(object){
    ll <- sum(sapply(object@GPs, function(gp) logLik(gp)))
    return(ll)
  }
)

#### summary ####
setMethod(
  f = "summary",
  signature = "GP_geomod",
  definition = function(object, ...){
    cat("Object of class ", class(object), "\n", sep = "")
    cat("Data points:", nrow(object@GPs[[1]]@data), "\n")
    Ntang <- ifelse(is.null(object@GPs[[1]]@tangents), 0,
                    nrow(object@GPs[[1]]@tangents))
    cat("Tangent points:", Ntang, "\n")
    cat("Number of classes:", length(object@GPs), "\n")
    cat("Total Log-likelihood:", logLik(object), "\n\n")
    lab <- names(object@GPs)
    for (i in seq_along(lab)){
      cat("* Class label:", lab[i], "\n")
      gp <- object@GPs[[i]]
      summary(gp)
    }
  }
)


#### Predict ####
#' @rdname Predict
setMethod(
  f = "Predict",
  signature = "GP_geomod",
  definition = function(object, target, to = "value", output.ind = T,
                        output.prob = T, use.unknown = T, Nsamp = 1e4){

    # setup
    ncat <- length(object@GPs)
    indmat <- varmat <- matrix(0, nrow(target), ncat)
    ydata <- GetData(target)
    class.names <- c(names(object@GPs), "Unknown")

    # prediction
    for (i in seq(ncat)){
      tmp <- Predict(object@GPs[[i]], target, to = "value", output.var = T)
      indmat[, i] <- tmp[["value"]]
      varstr <- ifelse(class(object@GPs[[i]]) == "SPGP",
                       "value.var_full", "value.var")
      varmat[, i] <- tmp[[varstr]] - object@GPs[[i]]@model@nugget
    }

    # probabilities
    randmat <- matrix(rnorm(Nsamp * ncat), Nsamp, ncat) # fix random numbers for smooth map
    prob <- t(sapply(seq(nrow(target)), function(i){
      tmpmat <- randmat * matrix(sqrt(varmat[i, ]), Nsamp, ncat, byrow = T) +
        matrix(indmat[i, ], Nsamp, ncat, byrow = T)
      tmpmat <- cbind(tmpmat, - rowSums(tmpmat))
      id <- apply(tmpmat, 1, which.max)
      as.numeric(table(factor(id, levels = 1:(ncat + 1)))) / Nsamp
    }))
    colnames(prob) <- paste0(to, "..", class.names, ".prob")

    # smoothing of probabilities to compensate for finite sample size
    prob <- t(apply(prob, 1, function(rw){
      Nzeros <- sum(rw == 0)
      rw <- rw * (1 - Nzeros / Nsamp)
      rw[rw == 0] <- 1 / Nsamp
      rw
    }))

    # skewed potential
    indmat_sk <- t(apply(log(prob), 1, function(x){
      xsort <- sort(x, decreasing = T)
      x - mean(xsort[1:2])
    }))
    inddf <- data.frame(indmat_sk)
    colnames(inddf) <- paste0(to, "..", class.names, ".ind")

    ## output
    if(!use.unknown){
      inddf <- inddf[, 1:ncat]
      prob <- prob[, 1:ncat]
      prob <- prob / matrix(rowSums(prob), nrow(prob), ncat)
    }

    # predicted class
    ydata[, to] <- apply(inddf, 1, function(rw) class.names[which.max(rw)])

    # indicators
    if(output.ind)
      ydata[, colnames(inddf)] <- inddf

    # probabilities
    if(output.prob){
      ydata[, colnames(prob)] <- prob
      ydata[, paste0(to, "..Entropy")] <-
        apply(prob, 1, function(rw) - sum(rw * log2(rw)))
    }

    # end
    target@data <- ydata
    return(target)
  }
)

#### Fit ####
#' @rdname Fit
# setMethod(
#   f = "Fit",
#   signature = "GP_geomod",
#   definition = function(object, contribution = T, nugget = T, maxrange = T,
#                         midrange = F, minrange = F, azimuth = F, dip = F, rake = F,
#                         power = F, ...){
#     lab <- names(object@GPs)
#     for(i in seq_along(object@GPs)){
#       # if(!(monitor == F)) cat("Fitting class", lab[i], "\n")
#       object@GPs[[i]] <- Fit(object@GPs[[i]], contribution = contribution,
#                              nugget = nugget, maxrange = maxrange,
#                              midrange = midrange, minrange = minrange,
#                              azimuth = azimuth, dip = dip, rake = rake,
#                              power = power, ...)
#     }
#
#     return(object)
#
#   }
# )
setMethod(
  f = "Fit",
  signature = "GP_geomod",
  definition = function(object, maxrange = T, midrange = F, minrange = F,
                        azimuth = F, dip = F, rake = F,
                        power = F, ...){

    # setup
    structures <- sapply(object@GPs[[1]]@model@structures, function(x) x@type)
    Nstruct <- length(structures)
    Ndata <- nrow(object@GPs[[1]]@data)
    data_box <- BoundingBox(object@GPs[[1]]@data)
    data_rbase <- sqrt(sum(data_box[1, ] - data_box[2, ])^2)
    GPs <- length(object@GPs)

    # optimization limits and starting point
    opt_min <- opt_max <- numeric(Nstruct * 8 + 1)
    xstart <- matrix(0, 1, Nstruct * 8 + 1)
    for (i in 1:Nstruct){
      # contribution
      opt_min[(i - 1) * 8 + 1] <- 0.1
      opt_max[(i - 1) * 8 + 1] <- 5
      xstart[(i - 1) * 8 + 1] <- object@GPs[[1]]@model@structures[[i]]@contribution

      # maxrange
      if (maxrange){
        opt_min[(i - 1) * 8 + 2] <- data_rbase / 1000
        opt_max[(i - 1) * 8 + 2] <- data_rbase * 5
      }
      else {
        opt_min[(i - 1) * 8 + 2] <- object@GPs[[1]]@model@structures[[i]]@maxrange
        opt_max[(i - 1) * 8 + 2] <- object@GPs[[1]]@model@structures[[i]]@maxrange
      }
      xstart[(i - 1) * 8 + 2] <- object@GPs[[1]]@model@structures[[i]]@maxrange

      # midrange (multiple of maxrange)
      if (midrange)
        opt_min[(i - 1) * 8 + 3] <- 0.01
      else
        opt_min[(i - 1) * 8 + 3] <- 1
      opt_max[(i - 1) * 8 + 3] <- 1
      xstart[(i - 1) * 8 + 3] <- object@GPs[[1]]@model@structures[[i]]@midrange /
        object@GPs[[1]]@model@structures[[i]]@maxrange

      # minrange(multiple of midrange)
      if (minrange)
        opt_min[(i - 1) * 8 + 4] <- 0.01
      else
        opt_min[(i - 1) * 8 + 4] <- 1
      opt_max[(i - 1) * 8 + 4] <- 1
      xstart[(i - 1) * 8 + 4] <- object@GPs[[1]]@model@structures[[i]]@minrange /
        object@GPs[[1]]@model@structures[[i]]@midrange

      # azimuth
      opt_min[(i - 1) * 8 + 5] <- 0
      if (azimuth)
        opt_max[(i - 1) * 8 + 5] <- 360
      else
        opt_max[(i - 1) * 8 + 5] <- 0
      xstart[(i - 1) * 8 + 5] <- object@GPs[[1]]@model@structures[[i]]@azimuth

      # dip
      opt_min[(i - 1) * 8 + 6] <- 0
      if (dip)
        opt_max[(i - 1) * 8 + 6] <- 90
      else
        opt_max[(i - 1) * 8 + 6] <- 0
      xstart[(i - 1) * 8 + 6] <- object@GPs[[1]]@model@structures[[i]]@dip

      # rake
      opt_min[(i - 1) * 8 + 7] <- 0
      if (rake)
        opt_max[(i - 1) * 8 + 7] <- 90
      else
        opt_max[(i - 1) * 8 + 7] <- 0
      xstart[(i - 1) * 8 + 7] <- object@GPs[[1]]@model@structures[[i]]@rake

      # power
      if (power){
        opt_min[(i - 1) * 8 + 8] <- 0.1
        opt_max[(i - 1) * 8 + 8] <- 3
      }
      else{
        opt_min[(i - 1) * 8 + 8] <- 1
        opt_max[(i - 1) * 8 + 8] <- 1
      }
      xstart[(i - 1) * 8 + 8] <- object@GPs[[1]]@model@structures[[i]]@power
    }

    # nugget
    opt_min[Nstruct * 8 + 1] <- 1e-6
    opt_max[Nstruct * 8 + 1] <- 2
    xstart[Nstruct * 8 + 1] <- 0.1 #mean(data_nugget)

    # conforming starting point to limits
    xstart[xstart < opt_min] <- opt_min[xstart < opt_min]
    xstart[xstart > opt_max] <- opt_max[xstart > opt_max]

    # fitness function
    makeGP <- function(x, finished = F){
      # covariance model
      m <- vector("list", Nstruct)
      for (i in 1:Nstruct){
        m[[i]] <- covarianceStructure3D(
          type = structures[i],
          contribution = x[( i - 1) * 8 + 1],
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
      model <- covarianceModel3D(x[Nstruct * 8 + 1], m)

      # GPs
      tmpgp <- object
      for(i in 1:GPs){
        tmpgp@GPs[[i]] <- GP(
          data = tmpgp@GPs[[i]]@data,
          model = model,
          value = "value",
          mean = tmpgp@GPs[[i]]@mean,
          trend = tmpgp@GPs[[i]]@trend,
          tangents = tmpgp@GPs[[i]]@tangents,
          force.interp = tmpgp@GPs[[i]]@data[["interpolate"]]
        )
      }
      # output
      tmpgp@params$nugget <- x[Nstruct * 8 + 1]
      if(finished)
        return(tmpgp)
      else
        return(logLik(tmpgp))
    }

    # optimization
    opt <- ga(
      type = "real-valued",
      fitness = function(x) makeGP(x, F),
      min = opt_min,
      max = opt_max,
      suggestions = xstart, ...
    )

    # update
    sol <- opt@solution
    return(makeGP(sol, T))
  }
)

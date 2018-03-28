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
#'
#' @seealso \code{\link{covarianceStructure3D-init}},
#' \code{\link{CovarianceMatrix}}
#'
#' @name covarianceModel3D-class
#' @export covarianceModel3D
covarianceModel3D <- setClass(
  "covarianceModel3D",
  slots = c(nugget = "numeric",
            nugget.dir = "numeric",
            structures = "list",
            total.var = "numeric"),
  validity = function(object){
    if (length(object@nugget) != 1)
      stop("nugget must be a scalar")
    if (length(object@nugget.dir) != 1)
      stop("nugget.dir must be a scalar")
    if (length(object@structures) > 1){
      if (!all(rapply(object@structures, class) == "covarianceStructure3D")){
        stop("multiple structures must all be of class 'covarianceStructure3D'")
      }
    }
    else if (!(class(object@structures[[1]]) %in%
               c("covarianceStructure3D",
                 "covarianceStructure3DNonStat"))){
      stop("single structure must be of class 'covarianceStructure3D'
           or 'covarianceStructure3DNonStat'")
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
#' @param structures An object of the \code{covarianceStructure3D} or
#' \code{covarianceStructure3DNonStat} class describing the desired covariance
#' structure.
#'
#' @details Stationary covariance accepts multiple structures, which can be
#' concatenated in a list. Non stationarity requires a single
#' \code{covarianceStructure3DNonStat} object.
#'
#' @name covarianceModel3D-init
#'
#' @seealso \code{\link{covarianceModel3D-class}},
#' \code{\link{covarianceStructure3D-class}}
covarianceModel3D <- function(nugget, structures, nugget.dir = 0){
  if (class(structures) != "list")
    structures <- list(structures)

  new("covarianceModel3D", nugget = nugget, nugget.dir = nugget.dir,
      structures = structures)
}

#### Getters ####
#' @rdname GetNonStatParams
setMethod(
  f = "GetNonStatParams",
  signature = "covarianceModel3D",
  definition = function(object, target){
    if (class(object@structures[[1]]) == "covarianceStructure3DNonStat")
      return(GetNonStatParams(object@structures[[1]], target))
    else
      stop("this model is stationary")
  }
)

#' @rdname GetPriorVariance
setMethod(
  f = "GetPriorVariance",
  signature = "covarianceModel3D",
  definition = function(object, target){
    if (class(object@structures[[1]]) == "covarianceStructure3DNonStat"){
      target <- Predict(object@structures[[1]]@contribution, target,
                        to = "contribution", output.var = F)
      return(exp(target[["contribution"]]))
    }
    else{
      # ajeitar para o linear
      total.var <- numeric(nrow(target))
      for (i in seq_along(object@structures)) {
        total.var <- total.var + GetPriorVariance(object@structures[[i]],
                                                  target)
      }
      return(total.var)
    }
  }
)

#### logLik ####
setMethod(
  f = "logLik",
  signature = "covarianceModel3D",
  definition = function(object){
    if (class(object@structures[[1]]) == "covarianceStructure3DNonStat")
      return(logLik(object@structures[[1]]))
    else
      return(0) # not really a log-likelihood
  }
)

#### show ####
setMethod(
  f = "show",
  signature = "covarianceModel3D",
  definition = function(object){

    if (class(object@structures[[1]]) == "covarianceStructure3DNonStat"){
      cat("Non stationary covariance model\n")
      cat("Type:", object@structures[[1]]@type, "\n")
      cat("Nugget:", object@nugget, "\n")
    }
    else{
      # statistics
      model_var <- sapply(object@structures, function(m) m@params$contribution)
      nug_rel <- 100 * object@nugget / (sum(model_var) + object@nugget)
      cov_rel <- 100 * model_var / (sum(model_var) + object@nugget)

      cat("Stationary covariance model\n")

      # covariance model
      cat("\nNugget: ", sprintf("%.4f", object@nugget),
          " (", sprintf("%02.2f", nug_rel),
          "%)\n\n", sep = "")
      for(i in seq_along(cov_rel)){
        txt <- paste0("Structure ", i, ": ",
                      sprintf("%.4f", object@structures[[i]]@params$contribution),
                      " (", sprintf("%02.2f", cov_rel[i]), "%)\n",
                      "Type: ", object@structures[[i]]@type, "\n")

        if (object@structures[[i]]@type != "bias")
          txt <- paste0(txt, "Range (max/med/min): ",
                        sprintf("%.2f", object@structures[[i]]@params$maxrange), "/",
                        sprintf("%.2f", object@structures[[i]]@params$midrange), "/",
                        sprintf("%.2f", object@structures[[i]]@params$minrange), "\n",
                        "Orientation (azimuth/dip/rake): ",
                        sprintf("%.2f", object@structures[[i]]@params$azimuth), "/",
                        sprintf("%.2f", object@structures[[i]]@params$dip), "/",
                        sprintf("%.2f", object@structures[[i]]@params$rake), "\n")

        if (object@structures[[i]]@type %in% c("linear", "nnet"))
          txt <- paste0(txt, "Center (X Y Z) = ",
                        paste0(sprintf("%.2f", object@structures[[i]]@params$center),
                               collapse = " "), "\n"
                        )

        if (object@structures[[i]]@type %in% c("cauchy", "nnet"))
          txt <- paste0(txt, "Power:",
                        sprintf("%.3f", object@structures[[i]]@params$power), "\n")

        cat(txt, "\n")
      }
    }

  }
)

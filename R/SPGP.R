#' @include GP.R
NULL

#### Sparse Pseudo-input Gaussian Process class ####
#' Sparse Gaussian Process
#'
#' Implementation of the Sparse Gaussian Process model for 3D spatial
#' interpolation. Extends the \code{GP} class.
#'
#' @slot data A \code{spatial3DDataFrame} object containing the necessary data.
#' @slot tangents A \code{directions3DDataFrame} object containing structural
#' geology data. Most likely generated with the \code{GetLineDirections()}
#' or \code{GetPlaneDirections()} method.
#' @slot model The covariance model. A \code{covarianceModel3D} object.
#' @slot mean The global mean. Irrelevant if a trend is used.
#' @slot trend The model's trend component. A formula in character format.
#' @slot beta The trend coefficients.
#' @slot likelihood The model's log-likelihood given the data.
#' @slot pre_comp A \code{list} containing pre-computed values to speed up
#' predictions.
#' @slot pseudo_inputs A \code{spatial3DDataFrame} object containing the coordinates
#' of the pseudo-inputs.
#' @slot pseudo_tangents A \code{directions3DDataFrame} object containing the coordinates
#' of the pseudo-inputs for directional data.
#' @slot variational A logical indicating if the model uses the variational
#' approach.
#'
#' @references
#' Snelson, E., Ghahramani, Z., 2006. Sparse Gaussian Processes using
#' Pseudo-inputs. Adv. Neural Inf. Process. Syst. 18 1257–1264.
#'
#' Titsias, M., 2009. Variational Learning of Inducing Variables in Sparse
#' Gaussian Processes. Aistats 5, 567–574.
#'
#' Bauer, M.S., van der Wilk, M., Rasmussen, C.E., 2016. Understanding
#' Probabilistic Sparse Gaussian Process Approximations. Adv. Neural Inf.
#' Process. Syst. 29.
#'
#' @seealso \code{\link{SPGP-init}}, \code{\link{GP-class}}
#' @name SPGP-class
#' @export SPGP
SPGP <- setClass(
  "SPGP",
  slots = c(pseudo_inputs = "spatial3DDataFrame",
            pseudo_tangents = "directions3DDataFrame",
            variational = "logical"
            ),
  contains = "GP"
)

#### initialization ####
#' Sparse Gaussian Process
#'
#' Implementation of the Sparse Gaussian Process model for 3D spatial
#' interpolation. Extends the \code{GP} class.
#'
#' @param data A \code{spatial3DDataFrame} object containing the data one
#' wishes to model.
#' @param tangents A \code{directions3DDataFrame} object containing structural
#' geology data. Most likely generated with the \code{GetLineDirections()}
#' method.
#' @param model The covariance model. A \code{covarianceStructure3D} or a
#' \code{list} containing one or more such objects.
#' @param value The column name of the variable to be modeled. It is assumed
#' the column does not contain missing values.
#' @param mean The global mean. Irrelevant if a trend is provided.
#' @param trend The model's trend component. A formula in character format.
#' @param pseudo_inputs The desired number of pseudo-inputs (whose coordinates
#' will be sampled from the data), a matrix or data frame with their
#' coordinates, or a 3D object.
#' @param pseudo_tangents The desired number of pseudo-structural data (whose
#' coordinates will be sampled from the data) or a \code{directions3DDataFrame}.
#' @param weights The importance of each data point in the model (a vector with
#' values between 0 and 1).
#' @param force.interp Indices of points that must be interpolated exactly.
#' @param reg.v Regularization to improve stability. A single value or a vector
#' with length matching the number of data points.
#' @param reg.t Regularization for structural data. A single value or a vector
#' with length matching the number of structural data.
#' @param variational Use the variational approach?
#'
#' @details This method builds a \code{SPGP} object with all the information
#' needed to make preditions at new data points.
#'
#' \code{trend} must be a character string with a formula as a function of
#' uppercase X, Y, and Z. The most common is the linear trend,
#' \code{"~ X + Y + Z"}. For ordinary kriging, use \code{"~1"}. If neither
#' \code{trend} nor \code{mean} are given, it is assumed that the global mean
#' is the mean of the data values.
#'
#' The SPGP works by compressing the information coming from all data into a
#' small number of pseudo-inputs. This way computational gains are obtained,
#' but the resulting model may pose difficulties for training.
#'
#' Given the sparse nature of the model, the points specified in \code{force.interp}
#' may still not be interpolated exactly. The effects of tangent data may also be
#' diminished.
#'
#' The variational model is cited in the literature as more stable and
#' less prone to overfitting. \code{variational = F} corresponds to the
#' Fully Independent Conditional (FIC) approach.
#'
#' @name SPGP-init
#'
#' @seealso \code{\link{SPGP-class}}, \code{\link{GP-class}}
#'
#' @references
#' Snelson, E., Ghahramani, Z., 2006. Sparse Gaussian Processes using
#' Pseudo-inputs. Adv. Neural Inf. Process. Syst. 18 1257–1264.
#'
#' Titsias, M., 2009. Variational Learning of Inducing Variables in Sparse
#' Gaussian Processes. Aistats 5, 567–574.
#'
#' Bauer, M.S., van der Wilk, M., Rasmussen, C.E., 2016. Understanding
#' Probabilistic Sparse Gaussian Process Approximations. Adv. Neural Inf.
#' Process. Syst. 29.
SPGP <- function(data, model, value,
                 mean = NULL, trend = NULL, pseudo_inputs,
                 weights = NULL,
                 force.interp = numeric(), reg.v = 1e-9,
                 tangents = NULL, reg.t = 1e-12,
                 pseudo_tangents = NULL, variational = T){

  # setup
  if (is.null(tangents))
    tangents <- as(tangents, "directions3DDataFrame")
  if (is.null(pseudo_tangents))
    pseudo_tangents <- as(pseudo_tangents, "directions3DDataFrame")
  if (is.null(trend) || (is.character(trend) & length(trend) == 0)){
    trend <- as.character(trend)
    beta <- matrix(0, 0, 0)
  }


  # pseudo-inputs
  if (class(pseudo_inputs) == "numeric" && nrow(data) >= pseudo_inputs) # number of inputs
    ps <- data[sample(nrow(data), pseudo_inputs), c("value")]
  else if(class(pseudo_inputs) %in% c("matrix", "data.frame")) # coordinates given
    ps <- points3DDataFrame(pseudo_inputs)
  else if (inherits(pseudo_inputs, "points3DDataFrame"))
    ps <- pseudo_inputs
  else
    stop("invalid format for pseudo-inputs")

  # pseudo-tangents
  if (!is.null(pseudo_tangents)){
    if(class(pseudo_tangents) == "numeric") # number of inputs
      pseudo_tangents <- tangents[sample(nrow(tangents), pseudo_tangents), ]
  }

  # covariance matrix of pseudo-inputs
  K_M <- Matrix(
    CovarianceMatrix(ps, ps, model) + diag(reg.v, nrow(ps), nrow(ps)),
    sparse = F
  )
  # pseudo-tangents
  if (!is.null(pseudo_tangents) && nrow(pseudo_tangents) > 0){
    K_M1 <- Matrix(CovarianceMatrix(ps, pseudo_tangents, model))
    K_M2 <- Matrix(CovarianceMatrix(pseudo_tangents, pseudo_tangents, model))
    K_M2 <- K_M2 + diag(reg.t, nrow(K_M2), ncol(K_M2))
    K_M <- rbind(
      cbind(K_M, K_M1),
      cbind(t(K_M1), K_M2)
    )
  }
  # symmetry and Cholesky decomposition
  K_M <- 0.5 * K_M + 0.5 * t(K_M)
  L <- t(chol(K_M))

  # if data is empty
  if (nrow(data) == 0){
    if(is.null(mean)) mean <- 0

    return(new("SPGP", data = data,
               tangents = as(tangents, "directions3DDataFrame"),
               model = model, mean = mean, variational = variational,
               likelihood = as.numeric(NA), trend = trend, beta = beta,
               pseudo_inputs = ps, pseudo_tangents = pseudo_tangents,
               pre_comp = list(LM = L, reg.v = reg.v, reg.t = reg.t)))
  }

  # value
  if(length(value) == 1 & class(value) == "character")
    data["value"] <- data[value]
  else
    data["value"] <- value
  yval <- data[["value"]]
  yval <- Matrix(yval, length(yval), 1)

  # weights
  if(is.null(weights)) weights <- rep(1, nrow(data))
  data["weights"] <- weights

  # interpolation
  int <- rep(F, nrow(data))
  int[force.interp] <- T
  data["interpolate"] <- int

  # mean
  if(is.null(mean)) mean <- mean(yval)

  # data
  data2 <- data[, c("value", "weights", "interpolate")]
  tangents <- as(tangents, "directions3DDataFrame")

  # covariance matrix of data and tangents
  K_NM <- Matrix(CovarianceMatrix(data, ps, model), sparse = F)
  if(!is.null(pseudo_tangents) && nrow(pseudo_tangents) > 0){
    K_NM <- cbind(K_NM,
                  Matrix(CovarianceMatrix(data, pseudo_tangents, model)))
  }
  if (!is.null(tangents) && nrow(tangents) > 0){

    K_tM <- t(Matrix(CovarianceMatrix(ps, tangents, model)))
    if(!is.null(pseudo_tangents) && nrow(pseudo_tangents) > 0){
      K_NM <- rbind(K_NM,
                    cbind(K_tM, Matrix(CovarianceMatrix(
                      tangents, pseudo_tangents, model
                    ))))
    }else{
      K_NM <- rbind(K_NM, K_tM)
    }
  }
  Ntang_ps <- nrow(pseudo_tangents)
  Ntang <- nrow(tangents)

  # approximated variance
  Q <- colSums(solve(L, t(K_NM))^2)

  # trend
  if(is.null(trend) | length(trend) == 0){
    H <- t(matrix(0, nrow(data) + Ntang, 0))
  }else{
    H <- t(TrendMatrix(data, trend))
    if(!is.null(tangents) && nrow(tangents) > 0){
      H <- cbind(H, t(TrendMatrix(tangents, trend)))
    }
    mean <- 0
  }
  Ntrend <- dim(H)[1]

  # pre-computations
  # mean
  yval <- yval - mean
  yval <- rbind(yval, Matrix(rep(0, Ntang), Ntang, 1))
  # delta
  delta <- model@total.var - Q
  if (Ntang_ps > 0)
    delta[-seq(nrow(data))] <- K_M2[1,1] - Q[-seq(nrow(data))]
  nugget_data <- rep(model@nugget, nrow(data))
  if (!is.null(weights)) nugget_data <- nugget_data / weights
  nugget_data[force.interp] <- reg.v
  nugget_data <- c(nugget_data, rep(reg.t, Ntang))
  if (variational)
    deltainv <- 1 /nugget_data
  else
    deltainv <- 1 / (delta + nugget_data)
  # matrices
  deltainv_mat <- Matrix(deltainv, nrow(data) + Ntang, nrow(ps) + Ntang_ps)
  dinvK_NM <- deltainv_mat * K_NM
  B <- K_M + t(K_NM) %*% dinvK_NM
  Breg <- c(rep(reg.v, nrow(data)), rep(reg.t, nrow(tangents)))
  B <- 0.5 * B + 0.5 * t(B) + diag(Breg, nrow(B), ncol(B)) # enforcing symmetry
  Binv <- solve(Matrix(B, sparse = F))
  K_Minv <- solve(Matrix(K_M))
  w_value <- Binv %*% t(K_NM) %*% (deltainv * yval)

  pre_comp <- list()
  pre_comp$w_value <- as.numeric(w_value)
  pre_comp$w_var <- as.matrix(K_Minv - Binv)
  pre_comp$B <- B
  pre_comp$LM <- L
  pre_comp$reg.v <- reg.v
  pre_comp$reg.t <- reg.t
  pre_comp$delta <- delta + nugget_data
  pre_comp$K_NM <- K_NM

  if(Ntrend > 0){
    F_ <- H %*% (matrix(deltainv, nrow(data) + Ntang, Ntrend) * t(H))
    Fi <- solve(F_)
    U <- H %*% dinvK_NM
    # G <- Fi - Fi %*% U %*% solve(-B + t(U) %*% Fi %*% U, t(U) %*% Fi)
    G <- F_ - U %*% Binv %*% t(U)
    J <- H %*% (deltainv * yval) - U %*% w_value
    beta <- solve(G, J)
    I <- diag(1, nrow(ps) + Ntang_ps, nrow(ps) + Ntang_ps)
    W_trend <- U %*% (Binv %*% (B-K_M) - I) %*% K_Minv

    beta <- as.matrix(beta)
    pre_comp$w_trend <- W_trend
    pre_comp$G <- G
  }


  # likelihood
  AL <- t(chol(model@nugget * B))
  gamma <- (delta + nugget_data) / nugget_data
  if (variational) gamma <- rep(1, length(gamma))
  ybar <- yval / sqrt(gamma)
  K_NMbar <- K_NM * Matrix(1 / sqrt(gamma), nrow(data) + Ntang, nrow(ps) + Ntang_ps)

  l1 <- - 0.5 * (sum(log(diag(AL))) - sum(log(diag(L))) + sum(gamma) +
                   (nrow(K_NM) - ncol(K_NM)) * log(model@nugget))
  l2 <- - 0.5 * (sum(ybar ^ 2) - sum(solve(AL, t(K_NMbar) %*% ybar) ^ 2)) / model@nugget


  # tmp <- dinvK_NM %*% Binv %*% (t(dinvK_NM) %*% yval)
  # l1 <- -0.5 * sum(yval * (deltainv * yval - tmp))
  # l2 <- -0.5 * (sum(delta + nugget_data) +
  #                 determinant(K_Minv)$modulus +
  #                 determinant(K_M + t(K_NM) %*%
  #                               dinvK_NM)$modulus)
  if (variational)
    l3 <- - 0.5 * sum(delta / nugget_data) # trace term
  else
    l3 <- 0
  lik <- l1 + l2 + l3 - 0.5 * length(yval) * log(2*pi)

  if(Ntrend > 0){
    l4 <- 0.5 * determinant(t(J) %*% beta)$modulus
    l5 <- 0.5 * determinant(G)$modulus
    lik <- lik + l4 - l5 + 0.5 * Ntrend * log(2*pi)
  }

  # end
  new("SPGP", data = data2,
      tangents = as(tangents, "directions3DDataFrame"),
      model = model, mean = mean, variational = variational,
      likelihood = as.numeric(lik), trend = as.character(trend), beta = beta,
      pseudo_inputs = ps, pseudo_tangents = pseudo_tangents,
      pre_comp = pre_comp)
}

#### show ####
setMethod(
  f = "show",
  signature = "SPGP",
  definition = function(object){
    # display
    cat("Object of class ", class(object), "\n", sep = "")
    if (object@variational)
      cat("Model type: Variational\n")
    else
      cat("Model type: FIC\n")
    cat("Data points:", nrow(object@data), "\n")
    cat("Pseudo-inputs:", nrow(object@pseudo_inputs), "\n")
    cat("Tangent points:", nrow(object@tangents), "\n")
    cat("Pseudo-tangents:", nrow(object@pseudo_tangents), "\n")
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
  signature = "SPGP",
  definition = function(object, target, to = "value", output.var = T){

    # trivial solution
    if(object@model@total.var == 0){
      target[, to] <- object@mean
      if(output.var){
        target[paste0(to, ".var_full")] <- 0
        target[paste0(to, ".var_cor")] <- 0
        target[paste0(to, ".quality")] <- 0
      }
      return(target)
    }

    # covariances
    ps <- object@pseudo_inputs
    Ktarget <- CovarianceMatrix(target, ps, object@model)
    if (!is.null(object@pseudo_tangents) && nrow(object@pseudo_tangents) > 0){
      Ktarget <- cbind(Ktarget,
                       CovarianceMatrix(target, object@pseudo_tangents,
                                        object@model))
    }

    # full variance
    tot_var <- object@model@total.var + object@model@nugget

    # correlated variance
    L <- object@pre_comp$LM
    Q_T <- colSums(solve(L, t(Ktarget))^2)

    # no data case
    if(nrow(object@data) == 0){
      target[, to] <- object@mean
      if(output.var){
        target[paste0(to, ".var_full")] <- tot_var
        target[paste0(to, ".var_cor")] <- Q_T
        target[paste0(to, ".quality")] <- Q_T / (tot_var - object@model@nugget)
      }
      return(target)
    }

    # trend
    if(length(object@trend) > 0){
      Htarget <- t(TrendMatrix(target, object@trend))
      w_tr <- object@pre_comp$w_trend
      G <- object@pre_comp$G
      beta <- object@beta
      R <- Htarget + w_tr %*% t(Ktarget)
    }

    # mean
    w_value <- object@pre_comp$w_value
    pred_mean <- apply(Ktarget, 1, function(rw) sum(rw * w_value))
    if(length(object@trend) > 0) pred_mean <- pred_mean + t(R) %*% beta
    target[to] <- pred_mean[1:nrow(target)] + object@mean

    # variance
    if(output.var){
      w_var <- object@pre_comp$w_var
      var_reduction <- rowSums((Ktarget %*% w_var) * Ktarget)[1:nrow(target)]

      # uncertainty due to trend
      var_trend <- numeric(nrow(target))
      if(length(object@trend) > 0){
        var_trend <- colSums(R * solve(G, R))
      }

      # full variance
      var_full <- tot_var - var_reduction + var_trend
      target[paste0(to, ".var_full")] <- var_full

      # correlated variance
      var_cor <- Q_T - var_reduction + var_trend
      target[paste0(to, ".var_cor")] <- var_cor

      # sparse approximation quality
      target[paste0(to, ".quality")] <- var_cor / (var_full - object@model@nugget)
    }

    # output
    return(target)
  }
)

#### Fit ####
#' @rdname Fit
setMethod(
  f = "Fit",
  signature = "SPGP",
  definition = function(object, contribution = F, maxrange = F,
                        midrange = F, minrange = F,
                        azimuth = F, dip = F, rake = F, pseudo_inputs = F,
                        pseudo_tangents = F,
                        power = F, nugget = F,
                        metric = c("logLik", "PLPD", "NRMSE"),
                        ...){
    require(GA)

    # setup
    structures <- sapply(object@model, function(x) x@type)
    Nstruct <- length(structures)
    Ndata <- nrow(object@data)
    Ntang <- nrow(object@tangents)
    data_var <- var(object@data[["value"]])
    data_nugget <- object@model@nugget
    Nps <- nrow(object@pseudo_inputs)
    Ntang_ps <- nrow(object@pseudo_tangents)
    metric <- metric[1]
    pseudo_inputs <- pseudo_inputs[1]

    # number of dimensions
    coords <- GetCoords(object@pseudo_inputs, "matrix")
    Ndims <- 3
    if (all(coords[, 3] == 0)) Ndims <- 2
    if (all(coords[, 2] == 0)) Ndims <- 1


    # bounding box
    data_box <- BoundingBox(object@data)
    ps_box <- BoundingBox(object@pseudo_inputs)
    data_box[1, ] <- apply(rbind(data_box[1, ], ps_box[1, ]), 2, min)
    data_box[2, ] <- apply(rbind(data_box[2, ], ps_box[2, ]), 2, max)
    if(Ntang > 0){
      tang_box <- BoundingBox(object@tangents)
      data_box[1, ] <- apply(rbind(data_box[1, ], tang_box[1, ]), 2, min)
      data_box[2, ] <- apply(rbind(data_box[2, ], tang_box[2, ]), 2, max)
    }
    data_rbase <- sqrt(sum(data_box[1,] - data_box[2,])^2)

    # continuous variables encoded with 10 bits (1024 values)
    bits <- c(rep(20, Nstruct * sum(contribution, maxrange, midrange,
                                       minrange, azimuth, dip, rake, power)),
              20 # nugget
              )#,
              # rep(10, Nps * 3)) # coordinates of pseudo-inputs

    # pseudo-inputs (subset of data)
    if (pseudo_inputs == "subset"){
      Nbits_ps <- ceiling(log2(Ndata / Nps + 0.001))
      bits <- c(bits, rep(Nbits_ps, Nps))
    }

    # pseudo-inputs (free to move)
    if (pseudo_inputs == "free")
      bits <- c(bits, rep(10, Nps * Ndims))

    # pseudo-tangents (subset of data)
    if (Ntang > 0 && pseudo_tangents){
      Nbits_tang <- ceiling(log2(Ntang / Ntang_ps + 0.001))
      bits <- c(bits, rep(Nbits_tang, Ntang_ps))
    }

    # fitness function
    makeGP <- function(x, finished = F){
      # decoding
      xdec <- .decodeString(x, bits)

      # covariance model
      m <- object@model@structures
      for(i in 1:Nstruct){
        tmp <- numeric(8)
        # contribution
        if (contribution){
          tmp[1] <- data_var * 10 ^ scales::rescale(xdec[1], from = c(0, 1048575),
                                                    to = c(-6, 2))
          xdec <- xdec[-1]
        }
        else
          tmp[1] <- m[[i]]@contribution
        # maxrange
        if (maxrange){
          tmp[2] <- data_rbase * 10 ^ scales::rescale(xdec[1], from = c(0, 1048575),
                                                      to = c(-3, 3))
          xdec <- xdec[-1]
        }
        else
          tmp[2] <- m[[i]]@maxrange
        # midrange
        if (midrange){
          tmp[3] <- scales::rescale(xdec[1], from = c(0, 1048575),
                                    to = c(0.01, 1)) * tmp[2]
          xdec <- xdec[-1]
        }
        else
          tmp[3] <- m[[i]]@midrange / m[[i]]@maxrange * tmp[2]
        # minrange
        if (minrange){
          tmp[4] <- scales::rescale(xdec[1], from = c(0, 1048575),
                                    to = c(0.01, 1)) * tmp[3]
          xdec <- xdec[-1]
        }
        else
          tmp[4] <- m[[i]]@minrange / m[[i]]@midrange * tmp[3]
        # azimuth
        if (azimuth){
          tmp[5] <- scales::rescale(xdec[1], from = c(0, 1048575),
                                    to = c(0, 360))
          xdec <- xdec[-1]
        }
        else
          tmp[5] <- m[[i]]@azimuth
        # dip
        if (dip){
          tmp[6] <- scales::rescale(xdec[1], from = c(0, 1048575),
                                    to = c(0, 90))
          xdec <- xdec[-1]
        }
        else
          tmp[6] <- m[[i]]@dip
        # rake
        if (rake){
          tmp[7] <- scales::rescale(xdec[1], from = c(0, 1048575),
                                    to = c(-90, 90))
          xdec <- xdec[-1]
        }
        else
          tmp[7] <- m[[i]]@rake
        # power
        if (power){
          tmp[8] <- scales::rescale(xdec[1], from = c(0, 1048575),
                                    to = c(0.1, 5))
          xdec <- xdec[-1]
        }
        else
          tmp[8] <- m[[i]]@power

        # building structure
        m[[i]] <- covarianceStructure3D(
          type = structures[i],
          contribution = tmp[1],
          maxrange = tmp[2],
          midrange = tmp[3],
          minrange = tmp[4],
          azimuth = tmp[5],
          dip = tmp[6],
          rake = tmp[7],
          power = tmp[8]
        )
      }

      # nugget
      if (nugget){
        nug <- data_var * 10 ^ scales::rescale(xdec[1], from = c(0, 1048575),
                                               to = c(-3, 2))
        xdec <- xdec[-1]
      }
      else
        nug <- object@model@nugget

      # pseudo-inputs (subset of data)
      if (pseudo_inputs == "subset"){
        ps_id <- .selectNofK(xdec[1:Nps] + 1, Ndata)
        xdec <- xdec[-(1:Nps)]
        ps <- GetCoords(object@data, "matrix")[ps_id, ]
      }
      # pseudo-inputs (free)
      else if (pseudo_inputs == "free"){
        ps <- matrix(xdec[1:(Ndims * Nps)], Nps, Ndims)
        ps <- cbind(ps, matrix(0, Nps, 3 - Ndims))
        for (i in 1:Ndims){
          ps[, i] <- scales::rescale(ps[, i], from = c(0, 1023),
                                     to = data_box[, i])
        }
        xdec <- xdec[-(1:(3 * Nps))]
      }
      else
        ps <- object@pseudo_inputs

      # pseudo_tangents (subset of data)
      if(Ntang_ps > 0 && pseudo_tangents){
        ps_tang_id <- .selectNofK(xdec[1:Ntang_ps] + 1, Ntang)
        ps_tang <- object@tangents[ps_tang_id, ]
      }
      else
        ps_tang <- object@pseudo_tangents

      # temporary GP
      tmpgp <- SPGP(
        data = object@data,
        model = covarianceModel3D(nug, m),
        value = "value",
        mean = object@mean,
        trend = object@trend,
        tangents = object@tangents,
        pseudo_inputs = ps,
        pseudo_tangents = ps_tang,
        reg.v = object@pre_comp$reg.v,
        reg.t = object@pre_comp$reg.t,
        force.interp = object@data[["interpolate"]],
        weights = object@data[["weights"]],
        variational = object@variational
      )
      # output
      if(finished)
        return(tmpgp)
      else{
        switch(metric,
               logLik = tmpgp@likelihood,
               PLPD = {
                 cv <- Xval(tmpgp)
                 cv$PLPD
               },
               NRMSE = {
                 cv <- Xval(tmpgp)
                 - cv$NRMSE
               }
               )
      }
    }

    # optimization
    opt <- ga(
      type = "binary",
      fitness = function(x) makeGP(x, F),
      nBits = sum(bits),
      ...
    )

    # update
    sol <- opt@solution[1, ]
    return(makeGP(sol, T))
  }
)

#### Simulate ####
#' @rdname Simulate
setMethod(
  f = "Simulate",
  signature = "SPGP",
  definition = function(object, target, Nsim, to = "value",
                        discount.noise = F, smooth = T, verbose = T,
                        randnum = NULL){

    # setup
    ps <- object@pseudo_inputs
    tot_var <- object@model@total.var
    if (nrow(object@data) == 0){
      B0 <- tcrossprod(object@pre_comp$LM)
      Bi0 <- solve(B0)
      w0 <- matrix(0, nrow(ps), 1)
    }
    else{
      B0 <- object@pre_comp$B
      Bi0 <- solve(B0)
      w0 <- B0 %*% as.matrix(object@pre_comp$w_value)
    }

    # random numbers
    if(is.null(randnum)){
      randmat <- matrix(rnorm(Nsim * nrow(target)), nrow(target), Nsim)
      pathmat <- matrix(0, nrow(target), Nsim)
      for(i in seq(Nsim)) pathmat[, i] <- sample(nrow(target))
    }

    else{
      randmat <- matrix(randnum, nrow(target), Nsim)
      pathmat <- matrix(seq(nrow(target)), nrow(target), Nsim)
    }

    # number formatting
    Ndigits <- length(unlist(strsplit(as.character(Nsim), NULL)))
    form <- paste0("%0", Ndigits, ".0f")

    # trivial solution
    if(object@model@total.var == 0){
      sim <- matrix(object@mean, nrow(target), Nsim)
      colnames(sim) <- paste0(to, ".sim_", sprintf(form, seq(Nsim)))
      target[, colnames(sim)] <- data.frame(as.matrix(sim))
      return(target)
    }

    # covariances
    K_TM <- Matrix(CovarianceMatrix(target, ps, object@model))
    L <- object@pre_comp$LM
    KMi <- as.matrix(solve(L %*% t(L)))
    Q_T <- colSums(solve(L, t(K_TM))^2)
    if(object@variational)
      d_T <- rep(object@pre_comp$reg.v, nrow(target))
    else{
      d_T <- tot_var - Q_T
    }


    # trend
    if(is.null(object@trend) | length(object@trend) == 0){
      H <- t(matrix(0, nrow(object@data), 0))
      yTR <- vTR <- matrix(0, nrow(target), 1)
    }else{
      H <- t(TrendMatrix(object@data, object@trend))
      H_T <- t(TrendMatrix(target, object@trend))

      R <- H_T + object@pre_comp$w_trend %*% t(K_TM)
      # contribution of trend to simulated values
      yTR <- t(R) %*% object@beta
      # contribution of trend to simulated variance
      vTR <- colSums(R * solve(object@pre_comp$G, R))
    }
    # if(discount.noise) vTR <- vTR - d_T

    # simulation
    for(i in seq(Nsim)){
      if(verbose) cat("\rSimulation", i, "of", Nsim, "...")
      useGPU = F
      if(useGPU){
        ysim <- .sparse_sim_gpu(path = pathmat[, i],
                                nugget = rep(object@model@nugget, nrow(target)),
                                w = as.numeric(w0),
                                Bi = as.matrix(Bi0),
                                KMi = KMi,
                                maxvar = rep(tot_var, nrow(target)),
                                K_TM = as.matrix(K_TM),
                                d_T = as.numeric(d_T),
                                yTR = as.numeric(yTR + object@mean),
                                vTR = as.numeric(vTR),
                                discount_noise = discount.noise,
                                Q_T = Q_T,
                                smooth = smooth,
                                randnum = randmat[, i])
      }
      else{
        ysim <- .sparse_sim(path = pathmat[, i],
                            nugget = rep(object@model@nugget, nrow(target)),
                            w_ = as.numeric(w0),
                            Bi_ = as.matrix(Bi0),
                            KMi_ = KMi,
                            maxvar = rep(tot_var, nrow(target)),
                            K_ = as.matrix(K_TM),
                            d_ = as.numeric(d_T),
                            yTR = as.numeric(yTR + object@mean),
                            vTR = as.numeric(vTR),
                            discount_noise = discount.noise,
                            Q = Q_T,
                            smooth = smooth,
                            randnum = randmat[, i])
      }
      target[paste0(to, ".sim_", sprintf(form, i))] <- ysim
    }
    if(verbose) cat("\n")
    return(target)
  }
)


#### Cross-validation ####
#' @rdname Xval
setMethod(
  f = "Xval",
  signature = "SPGP",
  definition = function(object){
    cv <- .SPGP_CV(
      as.numeric(object@pre_comp$B %*% as.matrix(object@pre_comp$w_value)),
      as.matrix(solve(object@pre_comp$B)),
      object@data[["value"]],
      as.matrix(object@pre_comp$K_NM),
      object@pre_comp$delta,
      rep(object@mean, nrow(object@data)),
      rep(0, nrow(object@data))
    )

    nug <- rep(object@nugget, nrow(object@data))
    nug[object@data[["interpolate"]]] <- object@pre_comp$reg.v
    d <- object@pre_comp$delta[1:nrow(object@data)]

    cv$RMSE <- sqrt(mean(((object@data[["value"]] - cv$mean) ^ 2)))
    cv$NRMSE <- sqrt(mean(((object@data[["value"]] - cv$mean) ^ 2) /
                            cv$var)) + sqrt(mean(d - nug))
    cv$LPD <- mean(dnorm(object@data[["value"]],
                         cv$mean, sqrt(cv$var), log = T))

    cv$PLPD <- cv$LPD - mean((d - nug) / nug)
    return(cv)
  }
)

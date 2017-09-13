#' @include GP.R
NULL

#### Sparse Pseudo-input Gaussian Process class ####
#' Sparse Gaussian Process
#'
#' A sparse model for working with low computational cost. Extends the
#' \code{GP} class.
#'
#' @slot data A \code{spatial3DDataFrame} object containing the necessary data.
#' @slot tangents A \code{directions3DDataFrame} object containing structural
#' geology data. Most likely generated with the \code{GetLineDirections()}
#' or \code{GetPlaneDirections()} method.
#' @slot model The covariance model. A \code{list} containing one or more
#' \code{covarianceStructure3D} objects.
#' @slot nugget The model's nugget effect or noise variance.
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
#' @references Snelson, E., Ghahramani, Z., 2006. Sparse Gaussian Processes
#' using Pseudo-inputs. Adv. Neural Inf. Process. Syst. 18 1257–1264.
#'
#' Titsias, M., 2009. Variational Learning of Inducing Variables in Sparse
#' Gaussian Processes. Aistats 5, 567–574.
#'
#' Bauer, M.S., van der Wilk, M., Rasmussen, C.E., 2016. Understanding
#' Probabilistic Sparse Gaussian Process Approximations. Adv. Neural Inf.
#' Process. Syst. 29.
#'
#' @seealso \code{\link{SPGP-init}}
#'
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
#' interpolation.
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
#' @param nugget The model's nugget effect or noise variance.
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
setMethod(
  f = "initialize",
  signature = "SPGP",
  definition = function(.Object, data, model, value, nugget,
                        mean = NULL, trend = NULL, pseudo_inputs,
                        weights = NULL,
                        force.interp = numeric(), reg.v = 1e-9,
                        tangents = NULL, reg.t = 1e-12,
                        pseudo_tangents = NULL, variational = T
                        ){
    require(Matrix)

    # covariance model
    if(length(model) == 1 & class(model) != "list")
      .Object@model <- list(model)
    else
      .Object@model <- model
    tot_var <- sum(sapply(.Object@model, function(md) md@contribution))

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
    .Object@data <- data[, c("value", "weights", "interpolate")]
    .Object@tangents <- as(tangents, "directions3DDataFrame")

    # pseudo-inputs
    if (class(pseudo_inputs) == "numeric") # number of inputs
      ps <- data[sample(nrow(data), pseudo_inputs), c("value")]
    else if(class(pseudo_inputs) %in% c("matrix", "data.frame")) # coordinates given
      ps <- points3DDataFrame(pseudo_inputs)
    else if (inherits(pseudo_inputs, "points3DDataFrame"))
      ps <- pseudo_inputs
    .Object@pseudo_inputs <- ps

    # pseudo-tangents
    if (!is.null(pseudo_tangents)){
      if(class(pseudo_tangents) == "numeric") # number of inputs
        pseudo_tangents <- tangents[sample(nrow(tangents), pseudo_tangents), ]
      .Object@pseudo_tangents <- pseudo_tangents
    }

    # covariances
    # pseudo-inputs
    K_M <- Matrix(
      CovarianceMatrix(ps, ps, model) + diag(reg.v, nrow(ps), nrow(ps))
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

    # data and tangents
    K_NM <- Matrix(CovarianceMatrix(data, ps, model))
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
    Ntang_ps <- nrow(.Object@pseudo_tangents)
    Ntang <- nrow(.Object@tangents)

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
    .Object@mean <- mean
    yval <- yval - mean
    yval <- rbind(yval, Matrix(rep(0, Ntang), Ntang, 1))
    # delta
    delta <- tot_var - Q
    if (Ntang_ps > 0)
      delta[-seq(nrow(data))] <- K_M2[1,1] - Q[-seq(nrow(data))]
    nugget_data <- rep(nugget, nrow(data))
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
    Breg <- c(rep(reg.v, nrow(data)), rep(reg.t, nrow(.Object@tangents)))
    B <- 0.5 * B + 0.5 * t(B) + diag(Breg, nrow(B), ncol(B)) # enforcing symmetry
    Binv <- solve(Matrix(B))
    K_Minv <- solve(Matrix(K_M))
    w_value <- Binv %*% t(K_NM) %*% (deltainv * yval)

    .Object@nugget <- nugget
    .Object@pre_comp$w_value <- as.numeric(w_value)
    .Object@pre_comp$w_var <- as.matrix(K_Minv - Binv)
    .Object@pre_comp$B <- B
    .Object@pre_comp$LM <- L
    .Object@pre_comp$reg.v <- reg.v
    .Object@pre_comp$reg.t <- reg.t
    .Object@pre_comp$delta <- delta + nugget_data
    .Object@pre_comp$K_NM <- K_NM
    .Object@variational <- variational

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

      .Object@beta <- as.matrix(beta)
      .Object@trend <- trend
      .Object@pre_comp$w_trend <- W_trend
      .Object@pre_comp$G <- G
    }


    # likelihood
    AL <- t(chol(nugget * B))
    gamma <- (delta + nugget_data) / nugget_data
    if (variational) gamma <- rep(1, length(gamma))
    ybar <- yval / sqrt(gamma)
    K_NMbar <- K_NM * Matrix(1 / sqrt(gamma), nrow(data) + Ntang, nrow(ps) + Ntang_ps)

    l1 <- - 0.5 * (sum(diag(AL) ^ 2) - sum(diag(L) ^ 2) + sum(gamma) +
                     (nrow(K_NM) - ncol(K_NM)) * log(nugget))
    l2 <- - 0.5 * (sum(ybar ^ 2) - sum(solve(AL, t(K_NMbar) %*% ybar) ^ 2)) / nugget


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

    .Object@likelihood <- as.numeric(lik)

    # end
    # validObject(.Object)
    return(.Object)
  }
)

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

    # covariances
    ps <- object@pseudo_inputs
    Ktarget <- CovarianceMatrix(target, ps, object@model)
    if (!is.null(object@pseudo_tangents) && nrow(object@pseudo_tangents) > 0){
      Ktarget <- cbind(Ktarget,
                       CovarianceMatrix(target, object@pseudo_tangents,
                                        object@model))
    }

    # full variance
    tot_var <- sum(sapply(object@model, function(m) m@contribution)) + object@nugget

    # correlated variance
    L <- object@pre_comp$LM
    Q_T <- colSums(solve(L, t(Ktarget))^2)

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
    pred_mean <- apply(Ktarget, 1, function(rw){
      sum(rw * w_value)
    })
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
    data_nugget <- object@nugget
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
      m <- object@model
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
        nug <- object@nugget

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
        model = m,
        value = "value",
        nugget = nug,
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

# setMethod(
#   f = "Fit",
#   signature = "SPGP",
#   definition = function(object, contribution = F, maxrange = F,
#                         midrange = F, minrange = F,
#                         azimuth = F, dip = F, rake = F, pseudo_inputs = F,
#                         pseudo_tangents = F,
#                         power = F, nugget = F, nugget_ps = F,
#                         nugget.fix = numeric(),
#                         monitor = F){
#     require(GA)
#
#     # setup
#     structures <- sapply(object@model, function(x) x@type)
#     Nstruct <- length(structures)
#     Ndata <- nrow(object@data)
#     Ntang <- nrow(object@tangents)
#     data_var <- var(object@data[["value"]])
#     data_nugget <- object@nugget
#     Nps <- nrow(object@pseudo_inputs)
#     Ntang_ps <- nrow(object@pseudo_tangents)
#
#     # bounding box
#     data_box <- BoundingBox(object@data)
#     if(Ntang > 0){
#       tang_box <- BoundingBox(object@tangents)
#       data_box[1, ] <- apply(rbind(data_box[1, ], tang_box[1, ]), 2, min)
#       data_box[2, ] <- apply(rbind(data_box[2, ], tang_box[2, ]), 2, max)
#     }
#     data_rbase <- sqrt(sum(data_box[1,] - data_box[2,])^2)
#
#     # optimization limits and starting point
#     opt_min <- opt_max <- numeric(Nstruct * 8 + 1 + Nps * 4 + Ntang_ps * 6)
#     xstart <- matrix(0, 1, Nstruct * 8 + 1 + Nps * 4 + Ntang_ps * 6)
#     for(i in 1:Nstruct){
#       # contribution
#       if(contribution){
#         opt_min[(i-1)*8+1] <- data_var / 1000
#         opt_max[(i-1)*8+1] <- data_var * 10
#       }else{
#         opt_min[(i-1)*8+1] <- object@model[[i]]@contribution
#         opt_max[(i-1)*8+1] <- object@model[[i]]@contribution
#       }
#       xstart[(i-1)*8+1] <- object@model[[i]]@contribution
#
#       # maxrange
#       if(maxrange){
#         opt_min[(i-1)*8+2] <- data_rbase / 100
#         opt_max[(i-1)*8+2] <- data_rbase * 2
#       }else{
#         opt_min[(i-1)*8+2] <- object@model[[i]]@maxrange
#         opt_max[(i-1)*8+2] <- object@model[[i]]@maxrange
#       }
#       xstart[(i-1)*8+2] <- object@model[[i]]@maxrange
#
#       # midrange (multiple of maxrange)
#       if(midrange)
#         opt_min[(i-1)*8+3] <- 0.01
#       else
#         opt_min[(i-1)*8+3] <- 1
#       opt_max[(i-1)*8+3] <- 1
#       xstart[(i-1)*8+3] <- object@model[[i]]@midrange /
#         object@model[[i]]@maxrange
#
#       # minrange(multiple of midrange)
#       if(minrange)
#         opt_min[(i-1)*8+4] <- 0.01
#       else
#         opt_min[(i-1)*8+4] <- 1
#       opt_max[(i-1)*8+4] <- 1
#       xstart[(i-1)*8+4] <- object@model[[i]]@minrange /
#         object@model[[i]]@midrange
#
#       # azimuth
#       opt_min[(i-1)*8+5] <- 0
#       if(azimuth)
#         opt_max[(i-1)*8+5] <- 360
#       else
#         opt_max[(i-1)*8+5] <- 0
#       xstart[(i-1)*8+5] <- object@model[[i]]@azimuth
#
#       # dip
#       opt_min[(i-1)*8+6] <- 0
#       if(dip)
#         opt_max[(i-1)*8+6] <- 90
#       else
#         opt_max[(i-1)*8+6] <- 0
#       xstart[(i-1)*8+6] <- object@model[[i]]@dip
#
#       # rake
#       opt_min[(i-1)*8+7] <- 0
#       if(rake)
#         opt_max[(i-1)*8+7] <- 90
#       else
#         opt_max[(i-1)*8+7] <- 0
#       xstart[(i-1)*8+7] <- object@model[[i]]@rake
#
#       # power
#       if(power){
#         opt_min[(i-1)*8+8] <- 0.1
#         opt_max[(i-1)*8+8] <- 3
#       }
#       else{
#         opt_min[(i-1)*8+8] <- 1
#         opt_max[(i-1)*8+8] <- 1
#       }
#       xstart[(i-1)*8+8] <- object@model[[i]]@power
#     }
#
#     # nugget
#     if(nugget){
#       opt_min[Nstruct * 8 + 1] <- data_var/1000
#       opt_max[Nstruct * 8 + 1] <- data_var*2
#     }
#     else{
#       opt_min[Nstruct * 8 + 1] <- 0 # not used
#       opt_max[Nstruct * 8 + 1] <- 0 # not used
#     }
#     xstart[Nstruct * 8 + 1] <- mean(data_nugget)
#
#     # pseudo_inputs
#     if(pseudo_inputs){
#       opt_min[(Nstruct * 8 + 2):(Nstruct * 8 + 1 + 3 * Nps)] <-
#         rep(data_box[1,], each = Nps)
#       opt_max[(Nstruct * 8 + 2):(Nstruct * 8 + 1 + 3 * Nps)] <-
#         rep(data_box[2,], each = Nps)
#     }
#     else{
#       opt_min[(Nstruct * 8 + 2):(Nstruct * 8 + 1 + 3 * Nps)] <-
#         as.numeric(object@pseudo_inputs)
#       opt_max[(Nstruct * 8 + 2):(Nstruct * 8 + 1 + 3 * Nps)] <-
#         as.numeric(object@pseudo_inputs)
#     }
#     xstart[(Nstruct * 8 + 2):(Nstruct * 8 + 1 + 3 * Nps)] <-
#       as.numeric(object@pseudo_inputs)
#
#     # pseudo_inputs (nugget)
#     pos <- Nstruct * 8 + 1 + 3 * Nps + 1
#     if(nugget_ps){
#       opt_min[pos:(pos + Nps - 1)] <- data_var/1000
#       opt_max[pos:(pos + Nps - 1)] <- data_var*2
#     }
#     else{
#       opt_min[pos:(pos + Nps - 1)] <- 0
#       opt_max[pos:(pos + Nps - 1)] <- 0
#     }
#     xstart[pos:(pos + Nps - 1)] <- object@nugget_ps
#
#     # pseudo_tangents
#     pos <- Nstruct * 8 + 1 + Nps * 4 + 1
#     if(pseudo_tangents){
#       opt_min[pos:(pos + 3 * Ntang_ps - 1)] <-
#         rep(data_box[1, ], each = Ntang_ps)
#       opt_max[pos:(pos + 3 * Ntang_ps - 1)] <-
#         rep(data_box[2,], each = Ntang_ps)
#       pos <- pos + Ntang_ps * 3
#       opt_min[pos:(pos + 3 * Ntang_ps - 1)] <- -1
#       opt_max[pos:(pos + 3 * Ntang_ps - 1)] <- 1
#     }
#     else{
#       opt_min[pos:(pos + 3 * Ntang_ps - 1)] <-
#         as.numeric(GetCoords(object@pseudo_tangents, "matrix"))
#       opt_max[pos:(pos + 3 * Ntang_ps - 1)] <-
#         as.numeric(GetCoords(object@pseudo_tangents, "matrix"))
#       pos <- pos + Ntang_ps * 3
#       opt_min[pos:(pos + 3 * Ntang_ps - 1)] <-
#         as.numeric(unlist(GetData(object@pseudo_tangents)))
#       opt_max[pos:(pos + 3 * Ntang_ps - 1)] <-
#         as.numeric(unlist(GetData(object@pseudo_tangents)))
#     }
#     pos <- Nstruct * 8 + 1 + Nps * 4 + 1
#     xstart[pos:(pos + 3 * Ntang_ps - 1)] <-
#       as.numeric(GetCoords(object@pseudo_tangents, "matrix"))
#     xstart[pos:(pos + 3 * Ntang_ps - 1)] <-
#       as.numeric(GetCoords(object@pseudo_tangents, "matrix"))
#     pos <- pos + Ntang_ps * 3
#     xstart[pos:(pos + 3 * Ntang_ps - 1)] <-
#       as.numeric(unlist(GetData(object@pseudo_tangents)))
#     xstart[pos:(pos + 3 * Ntang_ps - 1)] <-
#       as.numeric(unlist(GetData(object@pseudo_tangents)))
#
#     # conforming starting point to limits
#     xstart[xstart < opt_min] <- opt_min[xstart < opt_min]
#     xstart[xstart > opt_max] <- opt_max[xstart > opt_max]
#
#     # fitness function
#     makeGP <- function(x, finished = F){
#       # covariance model
#       m <- vector("list", Nstruct)
#       for(i in 1:Nstruct){
#         m[[i]] <- covarianceStructure3D(
#           type = structures[i],
#           contribution = x[(i-1)*8+1],
#           maxrange = x[(i-1)*8+2],
#           midrange = x[(i-1)*8+2] *
#             x[(i-1)*8+3],
#           minrange = x[(i-1)*8+2] *
#             x[(i-1)*8+3] *
#             x[(i-1)*8+4],
#           azimuth = x[(i-1)*8+5],
#           dip = x[(i-1)*8+6],
#           rake = x[(i-1)*8+7],
#           power = x[(i-1)*8+8]
#         )
#       }
#       # pseudo-inputs
#       ps <- matrix(
#         x[(Nstruct * 8 + 2):(Nstruct * 8 + 1 + 3 * Nps)],
#         Nps, 3)
#       pos <- Nstruct * 8 + 1 + 3 * Nps + 1
#       psnug <- x[pos:(pos + Nps - 1)]
#       # pseudo_tangents
#       if(Ntang_ps > 0){
#         pos <- Nstruct * 8 + 1 + Nps * 4 + 1
#         tang_coords <- matrix(x[pos:(pos + 3 * Ntang_ps - 1)], Ntang_ps, 3)
#         pos <- pos + Ntang_ps * 3
#         tang_df <- matrix(x[pos:(pos + 3 * Ntang_ps - 1)], Ntang_ps, 3)
#         tang_df <- apply(tang_df, 2, function(rw){
#           rw / sqrt(sum(rw ^ 2))
#         })
#         tang_df <- data.frame(tang_df)
#         colnames(tang_df) <- c("dX", "dY", "dZ")
#         ps_tang <- points3DDataFrame(tang_coords, tang_df)
#       }
#       else
#         ps_tang <- NULL
#
#       # temporary GP
#       if(nugget){
#         # fit a constant nugget model
#         tmpnug <- x[Nstruct * 8 + 1]
#       }
#       else{
#         # use values as given
#         tmpnug <- data_nugget
#       }
#       # points with fixed nugget
#       tmpnug[nugget.fix] <- data_nugget[nugget.fix]
#       # GP
#       tmpgp <- SPGP(
#         data = object@data,
#         model = m,
#         value = "value",
#         nugget = tmpnug,
#         mean = object@mean,
#         trend = object@trend,
#         tangents = object@tangents,
#         pseudo_inputs = ps,
#         nugget_ps = psnug,
#         pseudo_tangents = ps_tang,
#         reg.v = object@pre_comp$reg.v,
#         reg.t = object@pre_comp$reg.t
#       )
#       # output
#       if(finished)
#         return(tmpgp)
#       else
#         return(tmpgp@likelihood)
#     }
#
#     # optimization
#     opt <- ga(
#       type = "real-valued",
#       fitness = function(x) makeGP(x, F),
#       min = opt_min,
#       max = opt_max,
#       pmutation = 0.5,
#       popSize = 20,
#       run = 20,
#       monitor = monitor,
#       suggestions = xstart
#     )
#
#     # update
#     sol <- opt@solution
#     return(makeGP(sol, T))
#   }
# )

#### Simulate ####
#' @rdname Simulate
setMethod(
  f = "Simulate",
  signature = "SPGP",
  definition = function(object, target, Nsim, to = "value",
                        discount.noise = F, smooth = T, verbose = T){

    # setup
    B0 <- object@pre_comp$B
    Bi0 <- solve(B0)
    w0 <- B0 %*% as.matrix(object@pre_comp$w_value)
    ps <- object@pseudo_inputs
    tot_var <- sum(sapply(object@model, function(m) m@contribution))

    # number formatting
    Ndigits <- length(unlist(strsplit(as.character(Nsim), NULL)))
    form <- paste0("%0", Ndigits, ".0f")

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
      path <- sample(nrow(target))
      target[paste0(to, ".sim_", sprintf(form, i))] <-
        .sparse_sim(path = path,
                    nugget = object@nugget,
                    w_ = as.numeric(w0),
                    Bi_ = as.matrix(Bi0),
                    KMi_ = KMi,
                    maxvar = tot_var,
                    K_ = as.matrix(K_TM),
                    d_ = as.numeric(d_T),
                    yTR = as.numeric(yTR + object@mean),
                    vTR = as.numeric(vTR),
                    discount_noise = discount.noise,
                    Q = Q_T,
                    smooth = smooth)
    }
    cat("\n")
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

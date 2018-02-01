#### system stuff ####
.onUnload <- function (libpath) {
  library.dynam.unload("geomod3D", libpath)
}

#### GP auxiliary functions ####
.deriv_formula <- function(f, x){
  # get terms
  f2 <- attr(terms(f), "term.labels")
  if(length(f2) == 0) return(as.formula(paste0("~ I(0*",x,") - 1")))
  intercept <- attr(terms(f), "intercept") == 1
  # symbolic to arithmetic
  f2 <- gsub(":", "*", f2)
  f2 <- gsub(":", "*", f2)
  Is <- grep("I(.+)", f2)
  if(length(Is) > 0){
    f2[Is] <- gsub("[I(]", "", f2[Is])
    f2[Is] <- gsub(")$", "", f2[Is])
  }
  # derivative
  f2 <- sapply(f2, function(a) Deriv::Deriv(a, x))
  # arithmetic to symbolic
  zeros <- grep("^0$", f2)
  if(length(zeros) > 0){
    f2[zeros] <- paste0("0^",seq_along(zeros),"*",x)
  }
  f2 <- gsub("^1$", paste0(x,"^0"), f2)
  f2 <- paste0("I(", f2, ")", collapse = "+")
  if(intercept) f2 <- paste0("I(0*",x,") + ", f2)
  f2 <- paste("~", f2, "-1")
  f2 <- as.formula(f2)
}

# sparse sequential simulation using GPU
# .sparse_sim_gpu <- function(path, nugget, w, Bi, KMi, maxvar, K_TM,
#                             d_T, yTR, vTR, discount_noise, Q_T, smooth,
#                             randnum){
#
#   # conversion to GPU objects
#   w <- gpuMatrix(as.numeric(w), length(w), 1, type = "float")
#   Bi <- gpuMatrix(as.matrix(Bi), type = "float")
#   KMi <- gpuMatrix(as.matrix(KMi), type = "float")
#   # K_TM <- gpuMatrix(as.matrix(K_TM))
#   I <- gpuMatrix(diag(1, nrow(Bi), ncol(Bi)), type = "float")
#
#   vplus <- d_T + vTR
#
#   ysim <- matrix(rep(NA, length(path)))
#
#   # sequential simulation
#   for (i in seq_along(path)){
#     # k <- vclMatrix(K_TM[path[i], ], ncol(K_TM), 1)
#     k <- matrix(K_TM[path[i], ], ncol(K_TM), 1)
#
#     # prediction
#     mu <- (t(k) %*% (Bi %*% w))[1]
#     # v <- Q_T[path[i]] - (crossprod(k, ((KMi - Bi) %*% k)))[1] + vplus[path[i]]
#     v <- (t(k) %*% (Bi %*% k))[1] + vplus[path[i]]
#     if(v < 0) stop(paste("v =", v))
#
#     # simulation
#     ysim[path[i]] <- randnum[path[i]] * sqrt(v) + mu
#
#     # update
#     k <- k / sqrt(d_T[path[i]] + nugget[path[i]])
#     tmp <- 1 /(1 + crossprod(k, (Bi %*% k)))
#     Bi <- Bi %*% (I - k %*% crossprod(k, Bi)) * tmp[1]
#     k <- k / sqrt(d_T[path[i]] + nugget[path[i]])
#     w <- w + ysim[path[i]] * k
#   }
#
#   # smoothing
#   if (smooth){
#     # ysim <- gpuMatrix(ysim)
#     dnew <- matrix(1 / (maxvar - Q_T + 1e-9))
#     tmp <- matrix(rep(dnew[,], ncol(K_TM)), nrow(K_TM), ncol(K_TM))
#     tmp <- gpuMatrix(tmp, type = "float")
#     K_TM <- gpuMatrix(as.matrix(K_TM), type = "float")
#     Bnew <- gpuMatrix(solve(KMi)[])
#     Bnew <- Bnew + t(K_TM) %*% (tmp * K_TM)
#     Bnew <- Bnew + gpuMatrix(diag(1e-6, nrow(Bi), nrow(Bi)), type = "float") # regularization
#     ysim <- K_TM %*% gpuMatrix(solve(Bnew, crossprod(K_TM,  (ysim * dnew)))[]);
#     ysim <- ysim[];
#   }
#
#   # adding noise
#   if (!discount_noise)
#     ysim <- as.numeric(ysim) + rnorm(length(ysim), sd = sqrt(nugget))
#
#   # end
#   return(as.numeric(ysim) + yTR)
# }

#### drawing ####
.find_color_cont <- function(values, rng = range(values), col,
                             na.color = NA){
  valsc <- rescale(values, from = rng)
  cr <- colour_ramp(col, na.color = na.color)
  return(cr(valsc))
}

#### genetic algorithm ####
.decodeString <- function(string, bits){
  # Converts a binary string into a vector of integers. The 'bits' argument
  # is an integer vector containing the length of each segment.
  f <- unlist(sapply(seq_along(bits), function(i) rep(i, bits[i])))

  sapply(split(string, f), function(el){
    l <- as.logical(el)
    b <- binaryLogic::as.binary(l, logic = T)
    as.integer(b)
  })
}

.selectNofK <- function(string, K){
  # Selects a subset of  from a 1:K list. 'string' is an integer
  # vector containing the size of the "jumps" to be made. The function
  # will loop over the list as necessary.
  if (any(string == 0))
    stop("'string' argument must contain positive integers")

  if (length(string) == K)
    return(1:K)

  keep <- rep(F, K)
  pos <- 0
  for (i in seq_along(string)){
    pos <- pos + string[i]
    if (pos > K)
      pos <- 1
    if (pos %in% which(keep)){
      if (sum(which(!keep) > pos) > 0)
        pos <- which(!keep)[which(!keep) > pos][1]
      else
        pos <- max(which(!keep))
    }
    keep[pos] <- T
  }
  return(which(keep))
}

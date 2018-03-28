#' @include generics.R
#' @include utils.R
#' @include spatial3DDataFrame.R
#' @include points3DDataFrame.R
#' @include lines3DDataFrame.R
#' @include directions3DDataFrame.R
NULL

#### anisotropy ellipsoid ####
#' Anisotropy ellipsoid
#'
#' Builds a matrix that can be used to rotate and/or strech a set of coordinates
#' in order to model anisotropy.
#'
#' @param maxrange,midrange,minrange The  three semi-axes of the anisotropy
#' ellipsoid.
#' @param azimuth,dip,rake Orientation of the ellipsoid in geological
#' coordinates.
#' @param radians Are the angles in radians or degrees?
#'
#' @return A 3x3 matrix. Multiplication of a matrix of coordinates by this
#' matrix yields the transformed coordinates. Multiplication by its inverse
#' reverses the effect.
anisotropy3D <- function(maxrange, midrange = maxrange, minrange = midrange,
                         azimuth = 0, dip = 0, rake = 0, radians = F){
  # conversion to radians
  if(!radians){
    azimuth <- azimuth * pi/180
    dip <- dip * pi/180
    rake <- rake * pi/180
  }

  # conversion to mathematical coordinates
  dip <- -dip
  r <- c(midrange, maxrange, minrange)

  # rotation matrix
  Rx <- diag(1, 3, 3)
  Rx[c(1, 3), c(1, 3)] <- c(cos(rake), -sin(rake), sin(rake), cos(rake))
  Ry <- diag(1, 3, 3)
  Ry[2:3, 2:3] <- c(cos(dip), sin(dip), -sin(dip), cos(dip))
  Rz <- diag(1, 3, 3)
  Rz[1:2, 1:2] <- c(cos(azimuth), -sin(azimuth),
                   sin(azimuth), cos(azimuth))
  A <- diag(r, 3, 3)

  return(Rz %*% Ry %*% Rx %*% A)
}

#### covariance matrix ####
#' @rdname CovarianceMatrix
setMethod(
  f = "CovarianceMatrix",
  signature = c("lines3DDataFrame", "lines3DDataFrame", "covarianceModel3D"),
  definition = function(x, y = x, model, parts = 10){

    # line discretization
    xp <- Pointify(x, seq(0, 1, 1 / parts))
    yp <- Pointify(y, seq(0, 1, 1 / parts))

    # covariance matrix
    Kpoint <- CovarianceMatrix(xp, yp, model)
    K <- matrix(0, nrow(x), nrow(y))

    # averaging by pair of line segments
    for (i in seq(nrow(x))){
      for (j in seq(nrow(y))){
        idx <- ((i - 1) * parts + 1):(i * parts)
        idy <- ((j - 1) * parts + 1):(j * parts)
        K[i, j] <- mean(Kpoint[idx, idy])
      }
    }

    return(K)
  }
)

#' @rdname CovarianceMatrix
setMethod(
  f = "CovarianceMatrix",
  signature = c("lines3DDataFrame", "points3DDataFrame", "covarianceModel3D"),
  definition = function(x, y, model, parts = 10){

    # line discretization
    xp <- Pointify(x, seq(0, 1, 1 / parts))

    # covariance matrix
    Kpoint <- CovarianceMatrix(xp, y, model)
    K <- matrix(0, nrow(x), nrow(y))

    # averaging by line segment
    for (i in seq(nrow(x))){
      for (j in seq(nrow(y))){
        idx <- ((i - 1) * parts + 1):(i * parts)
        K[i, j] <- mean(Kpoint[idx, j])
      }
    }

    return(K)
  }
)

#' @rdname CovarianceMatrix
setMethod(
  f = "CovarianceMatrix",
  signature = c("points3DDataFrame", "lines3DDataFrame", "covarianceModel3D"),
  definition = function(x, y, model, parts = 10){
    CovarianceMatrix(y, x, model, parts)
  }
)

#' @rdname CovarianceMatrix
setMethod(
  f = "CovarianceMatrix",
  signature = c("points3DDataFrame","points3DDataFrame", "covarianceModel3D"),
  definition = function(x, y = x, model){
    # setup
    Nrows <- nrow(x)
    Ncols <- nrow(y)

    # covariance matrix
    u <- GetCoords(x, as = "matrix")
    v <- GetCoords(y, as = "matrix")
    K <- matrix(0, Nrows, Ncols)
    for (md in model@structures)
      K <- K + .CalcCovMat(md, u, v)
    return(K)

  }
)

#' @rdname CovarianceMatrix
setMethod(
  f = "CovarianceMatrix", # value/derivative covariance
  signature = c("points3DDataFrame", "directions3DDataFrame", "covarianceModel3D"),
  definition = function(x, y, model){
    # setup
    Ndata <- nrow(x)
    Ntang <- nrow(y)
    xcoords <- GetCoords(x, "matrix")
    tcoords <- GetCoords(y, "matrix")

    # dip and strike vectors
    vec <- y@directions

    # covariance matrix
    K <- matrix(0, Ndata, Ntang)
    for (md in model@structures)
      # K <- K + md@covd1(xcoords, tcoords, vec)
      K <- K +  .CalcCovMat_d1(md, xcoords, tcoords, vec)
    return(K)

  }
)

#' @rdname CovarianceMatrix
setMethod(
  f = "CovarianceMatrix", # derivative/value covariance
  signature = c("directions3DDataFrame", "points3DDataFrame", "covarianceModel3D"),
  definition = function(x, y, model){
    t(CovarianceMatrix(y, x, model))
  }
)

#' @rdname CovarianceMatrix
setMethod(
  f = "CovarianceMatrix", # derivative/derivative covariance
  signature = c("directions3DDataFrame", "directions3DDataFrame", "covarianceModel3D"),
  definition = function(x, y = x, model){
    # setup
    Ntang1 <- nrow(x)
    Ntang2 <- nrow(y)
    tcoords1 <- GetCoords(x, "matrix")
    tcoords2 <- GetCoords(y, "matrix")

    # tangent vectors
    dirvecs1 <- x@directions
    dirvecs2 <- y@directions

    # covariance matrix
    K <- matrix(0, Ntang1, Ntang2)
    for (md in model@structures) {
      # K <- K + md@covd2(tcoords1, tcoords2, dirvecs1, dirvecs2)
      K <- K + .CalcCovMat_d2(md, tcoords1, tcoords2, dirvecs1, dirvecs2)
    }
    return(K)

  }
)

#' @rdname CovarianceMatrix
setMethod(
  f = "CovarianceMatrix", # point/block covariance
  signature = c("points3DDataFrame", "blocks3DDataFrame", "covarianceModel3D"),
  definition = function(x, y, model){

    # setup
    Nrows <- nrow(x)
    Ncols <- nrow(y)
    yp <- Pointify(y)
    Ndisc <- prod(y@discretization)

    # covariance matrix
    u <- GetCoords(x, as = "matrix")
    v <- GetCoords(yp, as = "matrix")
    K <- matrix(0, Nrows, Ncols * Ndisc)
    for (md in model@structures)
      K <- K + .CalcCovMat(md, u, v)
      # K <- K + md@covfun(u, v)

    # averaging
    K <- t(rowsum(t(K), rep(seq(Ncols), each = Ndisc))) / Ndisc
    return(K)

  }
)

#' @rdname CovarianceMatrix
setMethod(
  f = "CovarianceMatrix", # block/point covariance
  signature = c("blocks3DDataFrame", "points3DDataFrame", "covarianceModel3D"),
  definition = function(x, y, model){
    t(CovarianceMatrix(y, x, model))
  }
)

#' @rdname CovarianceMatrix
setMethod(
  f = "CovarianceMatrix", # block/block covariance
  signature = c("blocks3DDataFrame", "blocks3DDataFrame", "covarianceModel3D"),
  definition = function(x, y = x, model){

    # setup
    Nrows <- nrow(x)
    Ncols <- nrow(y)
    xp <- Pointify(x)
    Ndiscx <- prod(x@discretization)
    yp <- Pointify(y)
    Ndiscy <- prod(y@discretization)

    # covariance matrix
    u <- GetCoords(xp, as = "matrix")
    v <- GetCoords(yp, as = "matrix")
    K <- matrix(0, Nrows * Ndiscx, Ncols * Ndiscy)
    for (md in model@structures)
      K <- K + .CalcCovMat(md, u, v)

    # averaging
    K <- t(rowsum(t(K), rep(seq(Ncols), each = Ndiscy))) / Ndiscy
    K <- rowsum(K, rep(seq(Nrows), each = Ndiscx)) / Ndiscx
    return(K)

  }
)

#' @rdname CovarianceMatrix
setMethod(
  f = "NuggetMatrix",
  signature = c("spatial3DDataFrame", "covarianceModel3D"),
  definition = function(x, model){

    K <- diag(model@nugget, nrow(x), nrow(x))
    return(K)
  }
)

#### trend matrix ####
#' @rdname TrendMatrix
setMethod(
  f = "TrendMatrix",
  signature = c("points3DDataFrame", "character"),
  definition = function(x, trend){
    trend <- as.formula(trend)
    TR <- model.matrix(trend, GetCoords(x, "data.frame"))
    return(TR)
  }
)

#' @rdname TrendMatrix
setMethod(
  f = "TrendMatrix",
  signature = c("directions3DDataFrame", "character"),
  definition = function(x, trend){
    # setup
    coords <- GetCoords(x, "data.frame")
    trend <- as.formula(trend)

    # trend gradient
    dTX <- model.matrix(.deriv_formula(trend, "X"), coords)
    dTY <- model.matrix(.deriv_formula(trend, "Y"), coords)
    dTZ <- model.matrix(.deriv_formula(trend, "Z"), coords)

    # directional vectors
    vec <- x@directions

    # directional derivatives
    dT <- dTX
    for (i in seq(dim(dTX)[2]))
      dT[, i] <- rowSums(cbind(dTX[, i], dTY[, i], dTZ[, i]) * vec)
    return(dT)
  }
)



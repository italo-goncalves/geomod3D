#' @include generics.R
#' @include utils.R
#' @include spatial3DDataFrame.R
#' @include points3DDataFrame.R
#' @include lines3DDataFrame.R
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
anisotropy3D <- function(maxrange, midrange, minrange = midrange,
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
  signature = c("lines3DDataFrame", "lines3DDataFrame"),
  definition = function(x, y, model, covariance = T, parts = 10){
    # covariance model
    if (length(model) == 1 & class(model) != "list"){
      model <- list(model)
    }
    if (!all(rapply(model, class) == "covarianceStructure3D")){
      stop("model must be of class 'covarianceStructure3D'")
    }

    # line discretization
    xp <- Pointify(x, seq(0, 1, 1 / parts))
    yp <- Pointify(y, seq(0, 1, 1 / parts))

    # covariance matrix
    Kpoint <- CovarianceMatrix(xp, yp, model, covariance)
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
  signature = c("lines3DDataFrame", "points3DDataFrame"),
  definition = function(x, y, model, covariance = T, parts = 10){
    # covariance model
    if (length(model) == 1 & class(model) != "list"){
      model <- list(model)
    }
    if (!all(rapply(model, class) == "covarianceStructure3D")){
      stop("model must be of class 'covarianceStructure3D'")
    }

    # line discretization
    xp <- Pointify(x, seq(0, 1, 1 / parts))

    # covariance matrix
    Kpoint <- CovarianceMatrix(xp, y, model, covariance)
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
  signature = c("points3DDataFrame", "lines3DDataFrame"),
  definition = function(x, y, model, covariance = T, parts = 10){
    CovarianceMatrix(y, x, model, covariance, parts)
  }
)

#' @rdname CovarianceMatrix
setMethod(
  f = "CovarianceMatrix",
  signature = c("points3DDataFrame","points3DDataFrame"),
  definition = function(x, y, model, covariance = T){
    # setup
    Nrows <- nrow(x)
    Ncols <- nrow(y)

    # covariance model
    if (length(model) == 1 & class(model) != "list"){
      model <- list(model)
    }
    if (!all(rapply(model, class) == "covarianceStructure3D")){
      stop("model must be of class 'covarianceStructure3D'")
    }

    # covariance matrix
    u <- GetCoords(x, as = "matrix")
    v <- GetCoords(y, as = "matrix")
    K <- matrix(0, Nrows, Ncols)
    if (covariance) for (md in model)
      K <- K + md@covfun(u, v)
    else for (md in model)
      K <- K + md@varfun(u, v)
    return(K)

  }
)

#' @rdname CovarianceMatrix
setMethod(
  f = "CovarianceMatrixD1", # value/derivative covariance
  signature = c("points3DDataFrame", "points3DDataFrame"),
  definition = function(x, tangents, model, covariance = T){
    # setup
    Ndata <- nrow(x)
    Ntang <- nrow(tangents)
    xcoords <- GetCoords(x, "matrix")
    tcoords <- GetCoords(tangents, "matrix")

    # covariance model
    if (length(model) == 1 & class(model) != "list"){
      model <- list(model)
    }
    if (!all(rapply(model, class) == "covarianceStructure3D")){
      stop("model must be of class 'covarianceStructure3D'")
    }

    # dip and strike vectors
    vec <- as.matrix(GetData(tangents[, c("dX","dY","dZ")]))

    # covariance matrix
    K <- matrix(0, Ndata, Ntang)
    if (covariance) for (md in model)
      K <- K + md@covd1(xcoords, tcoords, vec)
    else for (md in model)
      K <- K + md@vard1(xcoords, tcoords, vec)
    return(K)

  }
)

#' @rdname CovarianceMatrix
setMethod(
  f = "CovarianceMatrixD2", # derivative/derivative covariance
  signature = "points3DDataFrame",
  definition = function(tangents, model, covariance = T){
    # setup
    Ntang <- nrow(tangents)
    tcoords <- GetCoords(tangents, "matrix")

    # covariance model
    if (length(model) == 1 & class(model) != "list"){
      model <- list(model)
    }
    if (!all(rapply(model, class) == "covarianceStructure3D")){
      stop("model must be of class 'covarianceStructure3D'")
    }

    # tangent vectors
    dirvecs <- as.matrix(GetData(tangents[, c("dX","dY","dZ")]))

    # covariance matrix
    K <- matrix(0, Ntang, Ntang)
    if (covariance) for (md in model) {
      K <- K + md@covd2(tcoords, tcoords, dirvecs, dirvecs)
    }
    else for (md in model) {
      K <- K + md@vard2(tcoords, tcoords, dirvecs, dirvecs)
    }
    return(K)

  }
)

## trend matrix
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
  f = "TrendMatrixD1",
  signature = c("points3DDataFrame", "character"),
  definition = function(x, trend){
    # setup
    coords <- GetCoords(x, "data.frame")
    trend <- as.formula(trend)
    # trend gradient
    dTX <- model.matrix(.deriv_formula(trend, "X"), coords)
    dTY <- model.matrix(.deriv_formula(trend, "Y"), coords)
    dTZ <- model.matrix(.deriv_formula(trend, "Z"), coords)
    dTX <- rbind(dTX, dTX)
    dTY <- rbind(dTY, dTY)
    dTZ <- rbind(dTZ, dTZ)
    # dip and strike vectors
    vec <- as.matrix(GetData(x[, c("dX","dY","dZ")]))
    # directional derivatives
    dT <- dTX
    for (i in seq(dim(dTX)[2]))
      dT[, i] <- rowSums(cbind(dTX[, i], dTY[, i], dTZ[, i]) * vec)
    return(dT)
  }
)

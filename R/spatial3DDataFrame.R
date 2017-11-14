#' @include generics.R
NULL

#### spatial3DDataFrame class ####
#' Abstract class for 3D spatial objects
#'
#' @slot coords A \code{list} containing the coordinates for each data unit.
#' The contents of the \code{list} vary according to each object:
#' \itemize{
#'   \item For \code{points3DDataFrame} objects and its subclasses, a numeric
#'   vector of length 3;
#'   \item For \code{lines3DDataFrame} objects, a numeric vector of length 6,
#'   where the first 3 elements contain the coordinates of a line segment's
#'   start and the last 3 contain the coordinates of the segment's end.
#' }
#' @slot data A \code{data.frame} with each data unit's attributes.
#' @slot bbox A \code{matrix} containing the coordinates of two opposing
#' edges of the object's bounding box.
#' @details This class does not have a constructor method. It is used only to
#' define methods that are common to all of its subclasses.
#'
#' @seealso \code{\link{lines3DDataFrame-init}},
#' \code{\link{points3DDataFrame-init}}
spatial3DDataFrame <- setClass(
  "spatial3DDataFrame",
  slots = c(coords = "list",
            data = "data.frame",
            bbox = "matrix")
)

#### GetData ####
#' @rdname GetData
setMethod(
  f = "GetData",
  signature = "spatial3DDataFrame",
  definition = function(object){
    return(object@data)
  }
)

#### boundingBox ####
#' @rdname BoundingBox
setMethod(
  f = "BoundingBox",
  signature = "spatial3DDataFrame",
  definition = function(object){
    return(object@bbox)
  }
)

#### nrow, ncol ####
# data.frame analogs
setMethod("nrow", "spatial3DDataFrame",
          function(x){return(nrow(x@data))}
)
setMethod("ncol", "spatial3DDataFrame",
          function(x){return(ncol(x@data))}
)

#### as.data.frame ####
setMethod("as.data.frame", "spatial3DDataFrame",
          function(x){return(as(x, "data.frame"))}
)

#### show ####
setMethod(
  f = "show",
  signature = "spatial3DDataFrame",
  definition = function(object){
    # setup
    l <- min(10, nrow(object))
    if (l > 0){
      sp3df <- object[seq(l), ]
      coords <- GetCoords(sp3df,"data.frame")
      df <- GetData(sp3df)
      # display
      cat("Object of class ", class(object), "\n", sep = "")
      cat(nrow(object), " coordinates and ",
          ncol(object), " attributes\n\n", sep = "")
      cat("Bounding Box:\n")
      show(BoundingBox(object))
      cat("\nCoordinates:\n")
      show(coords)
      cat("\nAttributes:\n")
      show(df)
    }
    else
      cat("Empty", class(object))
  }
)

#### str ####
setMethod(
  f = "str",
  signature = "spatial3DDataFrame",
  definition = function(object, ...) str(GetData(object), ...)
)

#### [<- ####
setMethod(
  f = "[<-",
  signature = "spatial3DDataFrame",
  definition = function(x, i, j, value, drop){
    # checks
    if (missing(i)) i <- seq(nrow(x))
    if (class(i) == "character"){
      j <- i
      i <- seq(nrow(x))
    }
    if (class(j) != "character")
      stop(paste("Columns in", class(x),
                 "must be referenced by name"))
    if (class(value) %in% c("spatial3DDataFrame",
                            "points3DDataFrame",
                            "lines3DDataFrame")){
      value <- GetData(value)
    }
    if (!is.null(dim(value)) & class(value) != "data.frame")
      value <- as.matrix(value)
    # standardization as data frame
    value <- data.frame(value)
    colnames(value) <- j
    df <- GetData(x)
    # set on column at a time to avoid nasty bug
    for (k in j) df[i,k] <- value[,k]
    # output
    x@data <- df
    return(x)
  }
)

#### [[ ####
setMethod(
  f = "[[",
  signature = "spatial3DDataFrame",
  definition = function(x,i,exact){
    return(unlist(GetData(x)[,i]))
  }
)

#### GetNormals ####
#' @rdname GetNormals
setMethod(
  f = "GetNormals",
  signature = "spatial3DDataFrame",
  definition = function(object, dip = "Dip", strike = "Strike"){
    # setup
    df <- GetData(object)
    df[,".diprad"] <- -df[,dip] * pi/180
    df[,".strad"] <- (90 - df[,strike]) * pi/180
    # normal vector
    dipstr <- GetPlaneDirections(object, dip, strike)
    vec1 <- GetData(dipstr[1:nrow(object), ])
    vec2 <- GetData(dipstr[(nrow(object)+1):(2*nrow(object)), ])
    normalvec <- vec1[,c(2,3,1)] * vec2[,c(3,1,2)] -
      vec1[,c(3,1,2)] * vec2[,c(2,3,1)]
    # normalization
    normalvec <- t(apply(normalvec,1,function(v) v/sqrt(sum(v^2))))
    # result
    colnames(normalvec) <- c("nX", "nY", "nZ")
    normalvec <- as.data.frame(normalvec)
    return(directions3DDataFrame(coords = coords, directions = normalvec))
  }
)

#### GetPlaneDirections ####
#' @rdname GetPlaneDirections
setMethod(
  f = "GetPlaneDirections",
  signature = "spatial3DDataFrame",
  definition = function(object, dip = "Dip", strike = "Strike"){
    # setup
    df <- GetData(object)
    df[,".diprad"] <- -df[,dip] * pi/180
    df[,".strad"] <- (90 - df[,strike]) * pi/180
    # dip and strike
    dipvec <- matrix(c(
      cos(df[,".diprad"]) * cos(df[,".strad"] - pi/2),
      cos(df[,".diprad"]) * sin(df[,".strad"] - pi/2),
      sin(df[,".diprad"])
    ), nrow(object), 3)
    strvec <- matrix(c(
      cos(df[,".strad"]),
      sin(df[,".strad"]),
      rep(0, times = nrow(object))
    ), nrow(object), 3)
    # result
    vecs <- as.data.frame(rbind(dipvec, strvec))
    colnames(vecs) <- c("dX", "dY", "dZ")
    coords <- GetCoords(object)
    coords <- c(coords, coords)
    return(directions3DDataFrame(coords = coords, directions = vecs))
  }
)

#### GetLineDirections ####
#' @rdname GetLineDirections
setMethod(
  f = "GetLineDirections",
  signature = "spatial3DDataFrame",
  definition = function(object, dip = "Dip", azimuth = "Azimuth"){
    # setup
    df <- GetData(object)
    df[,".diprad"] <- -df[,dip] * pi/180
    df[,".azrad"] <- (90 - df[,azimuth]) * pi/180
    # vector directions
    dipvec <- matrix(c(
      cos(df[,".diprad"]) * cos(df[,".azrad"]),
      cos(df[,".diprad"]) * sin(df[,".azrad"]),
      sin(df[,".diprad"])
    ), nrow(object), 3)
    # result
    dipvec <- as.data.frame(dipvec)
    colnames(dipvec) <- c("dX", "dY", "dZ")
    coords <- GetCoords(object)
    return(directions3DDataFrame(coords = coords, directions = dipvec))
  }
)

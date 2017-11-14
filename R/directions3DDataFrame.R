#' @include points3DDataFrame.R

#### directions3DDataFrame class ####
#' 3D Directional data
#'
#' A formal representation of directional data. Extends the
#' \code{points3DDataFrame} class. It can be used to represent data on the
#' directional derivatives of a function.
#'
#' @slot coords A \code{list} containing the coordinates for each data unit.
#' @slot data A \code{data.frame} with each data unit's attributes.
#' @slot bbox A \code{matrix} containing the coordinates of two opposing
#' edges of the object's bounding box.
#' @slot directions A \code{matrix} containing coordinates of unit vectors
#' representing directions.
#'
#' @seealso \code{\link{directions3DDataFrame-init}},
#' \code{\link{spatial3DDataFrame-class}},
#' \code{\link{points3DDataFrame-init}}
#'
#' @export directions3DDataFrame
directions3DDataFrame <- setClass(
  "directions3DDataFrame",
  contains = "points3DDataFrame",
  slots = c(directions = "matrix")
)

#### initialization ####
#' 3D Directional data
#'
#' @param coords A list, matrix, or data frame containing the 3D coordinates
#' of the points.
#' @param df A data frame with the points' attributes
#' @param directions A matrix, or data frame containing the 3D coordinates
#' of the directional vectors.
#'
#' @details If \code{coords} or \code{directions} have less than 3
#' columns, the missing coordinates are given a value of 0. The directions are
#' rescaled to unit length.
#'
#' @seealso \code{\link{directions3DDataFrame-class}},
#' \code{\link{points3DDataFrame-class}}
#'
#' @name directions3DDataFrame
directions3DDataFrame <- function(coords, df, directions){

    # coordinates
    if (any(class(coords) %in% c("matrix", "Matrix", "data.frame", "tbl_df"))){
      coords <- as.matrix(coords)
      # enforcing 3D
      if (ncol(coords) > 3)
        stop("Invalid number of dimensions")
      coords <- cbind(coords, matrix(0, nrow(coords), 3 - ncol(coords)))
      # making list
      coords <- apply(coords, 1, function(x) list(x))
      coords <- lapply(coords, unlist)
    }else if (class(coords) == "list"){
      if (!(all(sapply(coords, length) == 3)))
        stop("Invalid number of dimensions")
    }else
      stop("Invalid format for coordinates")
    Ndata <- length(coords)

    # bounding box
    if (Ndata > 0){
      points_df <- data.frame(t(sapply(coords, function(z) z)))
      colnames(points_df) <- c("X", "Y", "Z")
      bbox <- as.matrix(rbind(
        apply(points_df, 2, min),
        apply(points_df, 2, max)))
      rownames(bbox) <- c("min", "max")
    }else
      bbox <- matrix(0, 2, 3)

    # directions
    if (any(class(directions) %in% c("matrix", "Matrix",
                                     "data.frame", "tbl_df"))){
      directions <- as.matrix(directions)
      # enforcing 3D
      if (ncol(directions) > 3)
        stop("Invalid number of dimensions")
      directions <- cbind(directions, matrix(0, nrow(directions),
                                             3 - ncol(directions)))
      # rescaling
      directions <- directions / matrix(sqrt(rowSums(directions ^ 2)),
                                        nrow(directions), 3)

      colnames(directions) <- c("dX", "dY", "dZ")
    }else
      stop("Invalid format for directions")

    # data
    if (missing(df)) df <- data.frame(.dummy = rep(NA, Ndata))

    # end
    new("directions3DDataFrame", coords = coords, data = df, bbox = bbox,
        directions = directions)
  }


#### setAs ####
setAs("NULL", "directions3DDataFrame", function(from, to)
  new(to, coords = list(), data = data.frame(), bbox = matrix(0, 0, 3),
      directions = matrix(0, 0, 3)))
setAs("directions3DDataFrame", "data.frame", function(from, to)
  cbind(GetCoords(from, "data.frame"), data.frame(from@directions),
        GetData(from)))

setMethod(
  f = "[",
  signature = "directions3DDataFrame",
  definition = function(x, i, j, drop){
    if (missing(i)) i <- seq(nrow(x))
    if (class(i) == "character"){
      j <- i
      i <- seq(nrow(x))
    }
    coords_list <- GetCoords(x)
    df <- GetData(x)
    coords_sub <- coords_list[i]
    df_sub <- df[i, j, drop = FALSE]
    dir_sub <- x@directions[i, j, drop = FALSE]
    return(new(class(x), coords = coords_sub, data = df_sub,
               directions = dir_sub))
  }
)

#### Pointify ####
#' @rdname Pointify
setMethod(
  f = "Pointify",
  signature = "directions3DDataFrame",
  definition = function(object){
    points3DDataFrame(object@coords, object@data)
  }
)

#### rbind, cbind equivalents ####
#' @rdname Bind
setMethod("Bind", c("directions3DDataFrame","directions3DDataFrame"),
          function(x, y){
            xp <- Pointify(x); yp <- Pointify(y)
            tmp <- Bind(xp, yp)
            directions3DDataFrame(tmp@coords, tmp@data,
                                  rbind(x@directions, y@directions))
          })

#### DrawDirections ####
#' @rdname DrawDirections
setMethod(
  f = "DrawDirections",
  signature = "directions3DDataFrame",
  definition = function(object, size, col = "yellow",
                        as = c("direction", "arrow")){
    # setup
    as = as[1]
    N <- nrow(object)
    coords <- GetCoords(object, "matrix")
    dirs <- object@directions
    if (length(col) < N)
      col <- rep(col, length.out = N)

    # drawing lines
    if (as == "direction"){
      for (i in seq(N)){
        cylcoords <- rbind(coords[i, ] + dirs[i, ] * size/2,
                           coords[i, ] - dirs[i, ] * size/2)
        shade3d(
          cylinder3d(cylcoords, radius = size/10,
                     sides = 32, closed = -2),
          col = col[i]
        )
      }
    }else if (as == "arrow"){
      for (i in seq(N))
        arrow3d(p0 = coords[i, ], p1 = coords[i, ] + dirs[i, ] * size,
                type = "rotation", col = col[i])
    }
  }
)


#' @include points3DDataFrame.R
NULL

#### grid3DDataFrame class ####
#' 3D regular point grid
#'
#' A 3D grid with equally spaced points. Extends the \code{points3DDataFrame}
#' class.
#'
#' @slot dims The number of grid points in each direction.
#'
#' @details Subsetting an object of this class will coerce the output to a
#' \code{points3DDataFrame} object.
#'
#' @seealso \code{\link{grid3DDataFrame-init}}
#'
#' @name grid3DDataFrame-class
#' @export grid3DDataFrame
grid3DDataFrame <- setClass(
  "grid3DDataFrame",
  slots = c(dims = "matrix"),
  contains = "points3DDataFrame",
  validity = function(object) {
    # callNextMethod(object)
  }
)

#### initialization ####
#' 3D regular point grid
#'
#' A 3D grid with equally spaced points. Extends the \code{points3DDataFrame}
#' class.
#'
#' @param gridx,gridy,gridz Grid coordinates in the three directions. The best
#' way to imput these is through the \code{seq()} function.
#' @param fields Column names to populate the grid. The columns are
#' automatically filled with \code{NA}s.
#'
#' @details The current implementation does not check if the points are
#' actually equally spaced in each direction. A future version will support
#' rotated grids.
#'
#' @name grid3DDataFrame-init
grid3DDataFrame <- function(gridx, gridy, gridz, fields = ".dummy"){
  # dimensions
  nx <- length(gridx)
  ny <- length(gridy)
  nz <- length(gridz)
  dims <- matrix(c(nx, ny, nz), 1, 3)
  colnames(dims) <- c("X","Y","Z")

  # coordinates
  coords <- c(
    rep(gridx, times = ny * nz),
    rep(gridy, each = nx, times = nz),
    rep(gridz, each = nx * ny)
  )
  coords <- matrix(coords, nx * ny * nz, 3)
  colnames(coords) <- c("X","Y","Z")

  # data
  nf <- length(fields)
  df <- data.frame(matrix(NA, nx * ny * nz, nf))
  colnames(df) <- fields

  # end
  p3df <- points3DDataFrame(coords, df)
  new("grid3DDataFrame", p3df, dims = dims)
}

#### show ####
setMethod(
  f = "show",
  signature = "grid3DDataFrame",
  definition = function(object){
    # setup
    l <- min(10, nrow(object))
    suppressWarnings(p3df <- object[seq(l), ])
    coords <- GetCoords(p3df, "data.frame")
    df <- GetData(p3df)
    # display
    cat("Object of class ", class(object), "\n", sep = "")
    cat(nrow(object), " coordinates and ",
        ncol(object), " attributes\n\n", sep = "")
    cat("Number of points:\n")
    d <- object@dims; rownames(d) <- ""; show(d)
    cat("\nBounding Box:\n")
    show(BoundingBox(object))
    cat("\nCoordinates:\n")
    show(coords)
    cat("\nAttributes:\n")
    show(df)
  }
)

#### [ ####
setMethod(
  f = "[",
  signature = "grid3DDataFrame",
  definition = function(x, i, j, drop){
    if (missing(i)) i <- seq(nrow(x))
    if (class(i) == "character"){
      j <- i
      i <- seq(nrow(x))
    }
    coords_list <- GetCoords(x)
    df <- GetData(x)
    coords_sub <- coords_list[i]
    df_sub <- df[i, j, drop=FALSE]
    return(points3DDataFrame(coords_sub, df_sub))
  }
)


#### Make3DArray ####
#' @rdname Make3DArray
setMethod(
  f = "Make3DArray",
  signature = "grid3DDataFrame",
  definition = function(object, value){
    df <- GetData(object)
    ar3d <- array(data = df[, value], dim = object@dims)
    coords <- GetCoords(object, "matrix")
    return(list(value = ar3d,
                x = sort(unique(coords[, 1])),
                y = sort(unique(coords[, 2])),
                z = sort(unique(coords[, 3])))
    )
  }
)

#### DrawSection ####
#' @rdname DrawSection
setMethod(
  f = "DrawSection",
  signature = "grid3DDataFrame",
  definition = function(object, by, values, col, col.default = "white",
                        x = NULL, y = NULL, z = NULL){

    # setup
    df <- GetData(object)
    if(!any(colnames(df) %in% by)){
      stop("Invalid attribute")
    }
    ar <- Make3DArray(object, by)
    d <- dim(ar$value)
    h <- c(ar$x[2] - ar$x[1],
           ar$y[2] - ar$y[1],
           ar$z[2] - ar$z[1]) * 0.5

    # pallette
    if (class(as.vector(ar$value)) == "numeric"){ # continuous variable
      objval <- object[[by]]
      colorsc <- .find_color_cont(objval, rng = range(values),
                                 col = col, na.color = col.default)
    }else{ # categorical variable
      # default color for categories not provided
      col <- c(col, col.default)
      ar$value.int <- rep(length(col) + 1, length(ar$value))
      for(i in seq_along(values)){
        ar$value.int[ar$value == values[i]] <- i
      }
      colorsc <- col[ar$value.int]
      dim(ar$value.int) <- d
      ar$value <- ar$value.int
    }
    dim(colorsc) <- d
    # plotting
    if (!is.null(x)){
      # nearest grid position
      id <- which.min(abs(x - ar$x))
      # color ordering
      section <- ar$value[id, , ]
      ord <- sort(section, index.return = T)
      colorsc <- colorsc[id, , ]
      colorsc <- unique(colorsc[ord$ix])
      # rescaling section
      v <- unique(sort(section))
      for (i in seq_along(v))
        section[section == v[i]] <- i
      # section
      show2d({
        par(mar=c(0, 0, 0, 0))
        image(section, col = colorsc)
        grid(d[2], d[3], col = "black", lty = "solid")
      },
      x = rep(ar$x[id], 4),
      y = c(ar$y[1] - h[2], ar$y[d[2]] + h[2],
            ar$y[d[2]] + h[2], ar$y[1] - h[2]),
      z = c(ar$z[1] - h[3], ar$z[1] - h[3],
            ar$z[d[3]] + h[3], ar$z[d[3]] + h[3]),
      width = 20 * d[2], height = 20 * d[3]
      )
      invisible(NULL)
    }else if (!is.null(y)){
      # nearest grid position
      id <- which.min(abs(y - ar$y))
      # color ordering
      section <- ar$value[, id, ]
      ord <- sort(section, index.return = T)
      colorsc <- colorsc[, id, ]
      colorsc <- unique(colorsc[ord$ix])
      # rescaling section
      v <- unique(sort(section))
      for (i in seq_along(v)){
        section[section == v[i]] <- i
      }
      # section
      show2d({
        par(mar=c(0, 0, 0, 0))
        image(section, col = colorsc)
        grid(d[1], d[3], col = "black", lty = "solid")
      },
      x = c(ar$x[1] - h[1], ar$x[d[1]] + h[1],
            ar$x[d[1]] + h[1], ar$x[1] - h[1]),
      y = rep(ar$y[id], 4),
      z = c(ar$z[1] - h[3], ar$z[1] - h[3],
            ar$z[d[3]] + h[3], ar$z[d[3]] + h[3]),
      width = 20 * d[1], height = 20 * d[3]
      )
      invisible(NULL)
    }else if (!is.null(z)){
      # nearest grid position
      id <- which.min(abs(z - ar$z))
      # color ordering
      section <- ar$value[, , id]
      ord <- sort(section, index.return = T)
      colorsc <- colorsc[, , id]
      colorsc <- unique(colorsc[ord$ix])
      # rescaling section
      v <- unique(sort(section))
      for(i in seq_along(v)){
        section[section == v[i]] <- i
      }
      # section
      show2d({
        par(mar=c(0, 0, 0, 0))
        image(section, col = colorsc)
        grid(d[1], d[2], col = "black", lty = "solid")
      },
      x = c(ar$x[1] - h[1], ar$x[d[1]] + h[1],
            ar$x[d[1]] + h[1], ar$x[1] - h[1]),
      y = c(ar$y[1] - h[2], ar$y[1] - h[2],
            ar$y[d[2]] + h[2], ar$y[d[2]] + h[2]),
      z = rep(ar$z[id], 4),
      width = 20 * d[1], height = 20 * d[2]
      )
      invisible(NULL)
    }else
      stop("One of x, y, or z must be provided")
  }
)

#### Pointify ####
#' @rdname Pointify
setMethod(
  f = "Pointify",
  signature = "grid3DDataFrame",
  definition = function(object){
    points3DDataFrame(GetCoords(object), GetData(object))
  }
)

#### SelectRegion ####
#' @rdname SelectRegion
setMethod(
  f = "SelectRegion",
  signature = "grid3DDataFrame",
  definition = function(object, xmin = -Inf, xmax = Inf,
                        ymin = -Inf, ymax = Inf, zmin = -Inf, zmax = Inf){
    # setup
    if (xmax <= xmin) stop("xmax must be grater than xmin")
    if (ymax <= ymin) stop("ymax must be grater than ymin")
    if (zmax <= zmin) stop("zmax must be grater than zmin")

    coords <- GetCoords(object, "matrix")

    # subsetting
    keep <- coords[, 1] >= xmin & coords[, 1] <= xmax &
      coords[, 2] >= ymin & coords[, 2] <= ymax &
      coords[, 3] >= zmin & coords[, 3] <= zmax

    # rebuilding
    df <- GetData(object)[keep, ]
    spacing <- apply(object@bbox, 2, diff) / object@dims
    coords_sub <- coords[keep, ]
    gr <- grid3DDataFrame(
      gridx = seq(min(coords_sub[, 1]), max(coords_sub[, 1]), spacing[1]),
      gridy = seq(min(coords_sub[, 2]), max(coords_sub[, 2]), spacing[2]),
      gridz = seq(min(coords_sub[, 3]), max(coords_sub[, 3]), spacing[3]),
      fields = colnames(df)
    )
    gr@data <- df

    return(gr)
  }
)

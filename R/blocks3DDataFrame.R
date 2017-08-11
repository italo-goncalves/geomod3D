#' @include grid3DDataFrame.R
NULL

#### blocks3DDataFrame class ####
#' 3D block model
#'
#' A 3D grid with adjacent blocks. Extends the \code{grid3DDataFrame}
#' class.
#'
#' @slot size The block size in the three dimensions.
#' @slot discretization The number of subdivisions for computation of
#' covariances.
#'
#' @details Subsetting an object of this class will coerce the output to a
#' \code{points3DDataFrame} object.
#'
#' @seealso \code{\link{blocks3DDataFrame-init}}
#'
#' @export blocks3DDataFrame
blocks3DDataFrame <- setClass(
  "blocks3DDataFrame",
  slots = c(size = "numeric",
            discretization = "numeric"),
  contains = "grid3DDataFrame"
)

#### initialization ####
#' 3D block model
#'
#' A 3D grid with adjacent blocks. Extends the \code{grid3DDataFrame}
#' class.
#'
#' @param corner A length 3 vector containing the smallest coordinate in each
#' direction.
#' @param number A length 3 vector containing the desired number of blocks in
#' each direction.
#' @param size A length 3 vector with the block dimensions in each direction.
#' @param fields Column names to populate the grid. The columns are
#' automatically filled with \code{NA}s.
#' @param discretization A length 3 vector indicating in how many parts to
#' "break" each block during the computations.
#'
#' @name blocks3DDataFrame-init
setMethod(
  f = "initialize",
  signature = "blocks3DDataFrame",
  definition = function(.Object, corner, number, size,
                        discretization = c(4, 4, 4),
                        fields = ".dummy"){
    # dimensions
    nx <- number[1]
    ny <- number[2]
    nz <- number[3]
    .Object@dims <- matrix(c(nx, ny, nz), 1, 3)
    colnames(.Object@dims) <- c("X","Y","Z")

    # block centers
    gridx <- seq(nx) * size[1] + corner[1] - size[1] / 2
    gridy <- seq(ny) * size[2] + corner[2] - size[2] / 2
    gridz <- seq(nz) * size[3] + corner[3] - size[3] / 2
    coords <- c(
      rep(gridx, times = ny * nz),
      rep(gridy, each = nx, times = nz),
      rep(gridz, each = nx * ny)
    )
    coords <- matrix(coords, nx * ny * nz, 3)
    colnames(coords) <- c("X","Y","Z")

    # bounding box
    .Object@bbox <- rbind(corner, corner + number * size)
    colnames(.Object@bbox) <- c("X", "Y", "Z")
    rownames(.Object@bbox) <- c("min", "max")

    # data
    nf <- length(fields)
    df <- data.frame(matrix(NA, nx * ny * nz, nf))
    colnames(df) <- fields

    # end
    p3df <- points3DDataFrame(coords,df)
    .Object@coords <- p3df@coords
    .Object@data <- p3df@data
    .Object@size <- size
    .Object@discretization <- discretization
    # validObject(.Object)
    return(.Object)
  }
)

#### show ####
setMethod(
  f = "show",
  signature = "blocks3DDataFrame",
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
    cat("Number of blocks:\n")
    d <- object@dims; rownames(d) <- ""; show(d)
    cat("Block size:", object@size, "\n")
    cat("Block discretization:", object@discretization, "\n")
    cat("\nBounding Box:\n")
    show(BoundingBox(object))
    cat("\nCoordinates:\n")
    show(coords)
    cat("\nAttributes:\n")
    show(df)
  }
)

#### Pointify ####
#' @rdname Pointify
setMethod(
  f = "Pointify",
  signature = "blocks3DDataFrame",
  definition = function(object, discretization = object@discretization){

    # discretized coordinates
    base_coords <- expand.grid(seq(discretization[1]),
                               seq(discretization[2]),
                               seq(discretization[3]))
    base_coords <- base_coords *
      matrix(object@size / discretization, nrow(base_coords), 3, byrow = T)
    base_coords <- base_coords -
      matrix(object@size / (1 * discretization), nrow(base_coords), 3, byrow = T)

    # full point grid
    center_coords <- GetCoords(object, "matrix")
    center_coords <- center_coords[rep(seq(nrow(object)),
                                       each = nrow(base_coords)), , drop = F]
    full_coords <- center_coords + base_coords[rep(seq(nrow(base_coords)),
                                                   nrow(object)), , drop = F]

    # data
    df <- GetData(object)[rep(seq(nrow(object)),
                              each = nrow(base_coords)), , drop = F]

    # end
    return(points3DDataFrame(full_coords, df))
  }
)

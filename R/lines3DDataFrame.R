#' @include spatial3DDataFrame.R

#### lines3DDataFrame class ####
#' Drillhole database
#'
#' Drillhole database represented as line segments. Extends the
#' \code{spatial3DdataFrame} class.
#'
#' @slot coords A \code{list} containing the coordinates for each data unit.
#' @slot data A \code{data.frame} with each data unit's attributes.
#' @slot bbox A \code{matrix} containing the coordinates of two opposing
#' edges of the object's bounding box.
#'
#' @details The \code{"HOLEID"} column is used by this object's methods. Do not
#' remove or rename it.
#'
#' @seealso \code{\link{lines3DDataFrame-init}},
#' \code{\link{spatial3DDataFrame-class}}
#'
#' @export lines3DDataFrame
lines3DDataFrame <- setClass(
  "lines3DDataFrame",
  contains = "spatial3DDataFrame",
  validity = function(object) {
    if (length(object@coords) !=
        nrow(object@data)) {
      stop("Number of coordinates does not match number of observations")
    }
    if (!all(rapply(object@coords,class) == "numeric")) {
      stop(
        "Invalid object in lines list. All objects must be of class 'numeric'"
      )
    }
  }
)

#### initialization ####
#' Drillhole database
#'
#' Drillhole database represented as line segments.
#'
#' @param lines_list A \code{list} containing the line segments' coordinates.
#' @param df A \code{data.frame} with the line segments' attributes.
#' @param collar,assay,survey Data frames with drillhole data.
#' @param holeid,from,to,X,Y,Z Column names from which to draw the parameters.
#'
#' @details There are two ways to build a \code{lines3DDataFrame} object:
#' either providing the data in the format that is seen in the final object
#' (less common) or by providing collar and assay data (more common).
#'
#' In the first mode, both \code{lines_list} and \code{df} must be given.
#' \code{lines_list} is a \code{list} where each entry is a numeric vector
#' of length 6 containing the X, Y, and Z coordinates for both ends of a
#' line segment. \code{df} is a data frame containing the attributes for
#' each segment.
#'
#' In the second mode, both \code{collar} and \code{assay} must be given.
#' \code{collar} contains the X, Y, and Z coordinates for each drillhole's
#' collar. \code{assay} contains the lithological description, chemical
#' analysis, etc. in intervals delimimted by the \code{from} and \code{to}
#' columns. \code{survey} contains the dipmeter measurements for each interval.
#' If it is missing, the holes are assumed vertical. The \code{holeid} column
#' must be present in all three data frames.
#'
#' @seealso \code{\link{spatial3DDataFrame-class}},
#' \code{\link{lines3DDataFrame-class}}
#'
#' @name lines3DDataFrame-init
setMethod(
  f = "initialize",
  signature = "lines3DDataFrame",
  definition = function(.Object, lines_list = NULL, df = NULL,
                        collar = NULL, assay = NULL, survey = NULL,
                        holeid = "HOLEID", from = "FROM", to = "TO",
                        X = "X", Y = "Y", Z = "Z"){

    ## building from assay data
    if (!is.null(collar)){
      # setup
      collar_names <- colnames(collar)
      collar_names <- collar_names[!(collar_names %in% c(X,Y,Z))]
      assay_names <- colnames(assay)
      assay_names <- assay_names[!(assay_names %in% c(holeid,from,to))]

      # processing data
      df <- merge(collar, assay, by=holeid)
      lines_list <- vector("list", nrow(df))
      if(is.null(survey)){
        for(i in seq(nrow(df))){
          lines_list[[i]] <- c(
            df[i, c(X,Y,Z)] - c(0,0,df[i,from]),
            df[i, c(X,Y,Z)] - c(0,0,df[i,to])
          )
        }
        lines_list <- lapply(lines_list, unlist)
      }else{
        # fazer survey aqui
        stop("Survey data not supported yet")
      }
      df <- df[, c(collar_names, assay_names)]
    }

    ## building directly from list and data.frame
    .Object@coords <- lines_list
    .Object@data <- df
    # bounding box
    lines_df <- GetCoords(.Object,"data.frame")
    names(lines_df) <- rep(c("X","Y","Z"),2)
    lines_df <- rbind(lines_df[,1:3],lines_df[,4:6])
    bbox <- as.matrix(rbind(
      apply(lines_df,2,min),
      apply(lines_df,2,max)))
    rownames(bbox) <- c("min","max")
    .Object@bbox <- bbox
    validObject(.Object)
    return(.Object)
  }
)

#### getCoords ####
#' @rdname GetCoords
setMethod(
  f = "GetCoords",
  signature = "lines3DDataFrame",
  definition = function(object, as = c("list", "data.frame", "matrix")){
    as <- as[1]
    if (as == "list"){
      return(object@coords)
    }else if (as == "data.frame"){
      lines_list <- object@coords
      df <- t(sapply(lines_list, function(x) x))
      df <- as.data.frame(df)
      colnames(df) <- c("X.1","Y.1","Z.1","X.2","Y.2","Z.2")
      return(df)
    }else if (as == "matrix"){
      lines_list <- object@coords
      df <- t(sapply(lines_list, function(x) x))
      colnames(df) <- c("X.1","Y.1","Z.1","X.2","Y.2","Z.2")
      return(df)
    }else{
      stop("Invalid value for 'as' argument")
    }
  }
)


#### getLength ####
#' @rdname GetLength
setMethod(
  f = "GetLength",
  signature = "lines3DDataFrame",
  definition = function(object){
    lines_list <- GetCoords(object)
    l <- sapply(lines_list, function(x){
      sqrt(sum((x[1:3] - x[4:6])^2))
    })
    return(l)
  }
)


#### show ####
setMethod(
  f = "show",
  signature = "lines3DDataFrame",
  definition = function(object){
    # setup
    l <- min(10, nrow(object))
    l3df <- object[seq(l), ]
    coords <- GetCoords(l3df,"data.frame")
    df <- GetData(l3df)
    # display
    cat("Object of class ", class(object), "\n", sep = "")
    cat(nrow(object), " line segments and ",
        ncol(object), " attributes\n\n", sep = "")
    cat("Bounding box:\n")
    show(BoundingBox(object))
    cat("\nLine segments:\n")
    show(coords)
    cat("\nAttributes:\n")
    show(df)
  }
)

#### Pointify ####
#' @rdname Pointify
setMethod(
  f = "Pointify",
  signature = "lines3DDataFrame",
  definition = function(object, locations = c(0.05, 0.5, 0.95), distance = F){
    if (min(locations) < 0 || max(locations) > 1)
      stop("locations must contain values between 0 and 1, inclusive")

    # function to break line into points
    break_line <- function(line, locations){
      lstart <- line[1:3]
      lend <- line[4:6]
      ldif <- lend - lstart
      points_list <- vector("list", length(locations))
      for(i in seq_along(locations)){
        points_list[[i]] <- lstart + locations[i] * ldif
      }
      return(points_list)
    }

    # conversion to points3DDataFrame
    r <- length(locations)
    n <- nrow(object)
    d <- numeric(n * r)
    points_list <- vector("list", n * r)
    lines_list <- GetCoords(object)
    l <- GetLength(object)
    df <- GetData(object)
    for (i in seq_len(n)){
      points_list[(r*(i-1)+1):(i*r)] <-
        break_line(lines_list[[i]], locations)
      d[(r*(i-1)+1):(i*r)] <- apply(rbind(
        l[i] * locations,
        l[i] * (1-locations)
      ), 2, min)
    }
    df <- df[rep(seq(n), each = r),,drop=F]
    if (distance) df[,".dist"] <- d
    return(points3DDataFrame(points_list,df))
  }
)

#### DrawDrillholes ####
#' @rdname DrawDrillholes
setMethod(
  f = "DrawDrillholes",
  signature = "lines3DDataFrame",
  definition = function(object, by, values, col, size, col.default = "white"){

    # setup
    object <- MergeSegments(object, by)
    df <- GetData(object)
    lines_list <- GetCoords(object)
    if (!any(colnames(df) %in% by)){
      stop("Invalid attribute")
    }
    N <- nrow(object)
    if (length(size) < N) size <- rep(size, length.out = N)
    objval <- df[, by]

    # pallette
    if(class(objval) == "numeric"){ # continuous variable
      colorsc <- .find_color_cont(objval,
                                 rng = range(values),
                                 col = col, na.value = col.default)
    }else{ # categorical variable
      names(col) <- values
      colorsc <- col[objval]
    }
    # plotting
    for(i in seq(N)){
      shade3d(
        cylinder3d(rbind(lines_list[[i]][1:3], lines_list[[i]][4:6]),
                   sides = 16, radius = size[i]/2,
                   closed = -2),
        col = colorsc[i])
    }
  }
)

#### DrawHoleID ####
#' @rdname DrawHoleID
setMethod(
  f = "DrawHoleID",
  signature = "lines3DDataFrame",
  definition = function(object, cex = 1){
    # setup
    coords <- GetCoords(object, "data.frame")
    df <- cbind(coords, GetData(object))
    holeids <- unique(df$HOLEID)
    N <- length(holeids)
    coords2 <- as.data.frame(matrix(NA, N, 3))

    # highest Z
    for(i in seq(N)){
      tmp <- df[df$HOLEID == holeids[i], ]
      pos <- which(tmp$Z.1 == max(tmp$Z.1))
      coords2[i, ] <- tmp[pos, 1:3]
    }

    # plotting
    text3d(coords2, texts = holeids, cex = cex, adj = c(0.5, 0))

  }
)

#### MergeSegments ####
#' @rdname MergeSegments
setMethod(
  f = "MergeSegments",
  signature = "lines3DDataFrame",
  definition = function(object, by,
                        keep = colnames(GetData(object))){

    # setup
    line_lengths <- GetLength(object)
    line_coords <- GetCoords(object, "matrix")
    lines_list <- GetCoords(object)
    df <- GetData(object)
    if (by == "HOLEID"){
      by <- ".HOLEID"
      df[,by] <- df[,"HOLEID"]
    }
    if (any(keep %in% c("HOLEID", by))){
      keep <- keep[-which(keep %in% c("HOLEID", by))]
    }
    df <- df[,c("HOLEID", by, keep)]
    Nlines <- nrow(line_coords)
    directions <- (line_coords[, 1:3] - line_coords[, 4:6]) / line_lengths

    # finding mergeable segments
    # condition 1 - segments sharing a point
    coord_from <- line_coords[2:Nlines, 1:3]
    coord_to <- line_coords[1:(Nlines-1), 4:6]
    start_end <- apply(coord_from - coord_to, 1, function(x) all(x==0))
    # condition 2 - parallelism
    dir_from <- directions[2:Nlines, ]
    dir_to <- directions[1:(Nlines-1), ]
    parallel <- apply(dir_from - dir_to, 1, function(x) all(x==0))
    # condition 3 - same value in "by" column
    val_to <- df[1:(Nlines-1), by]
    val_from <- df[2:Nlines, by]
    same_value <- val_from == val_to
    # final vector
    mergeable <- start_end & parallel & same_value

    # building merged object
    # find contiguous mergeable segments
    merge_ids <- split(seq(Nlines), diffinv(!mergeable))
    # number of rows in new object
    Nlines2 <- length(merge_ids)
    # new data frame
    df2 <- df[, -which(colnames(df) %in% c("HOLEID", by)), drop=F]
    cols_to_pass <- colnames(df2)
    df2[, ".line"] <- seq(nrow(df2)) # unique ID to avoid spread error
    for(i in cols_to_pass){
      if(class(df2[, i]) != "numeric"){
        df2[, i] <- as.character(df2[, i])
        # special characters in column names may give error here
        df2[, i] <- paste0(i, "..", make.names(df2[, i]))
        df2[, ".val"] <- 1
        df2 <- tidyr::spread_(data = df2, key_col = i, value_col = ".val",
                              fill = 0)
        df2 <- dplyr::arrange_(df2, ".line")
      }
    }
    df2 <- df2[, -which(colnames(df2) == ".line")]
    df2 <- cbind(df[, c("HOLEID", by), drop=F], df2)
    # averaging values
    lines_list2 <- vector("list", Nlines2)
    df3 <- data.frame(matrix(NA, Nlines2, ncol(df2)))
    colnames(df3) <- colnames(df2)
    for (i in seq_along(merge_ids)){
      l <- length(merge_ids[[i]])
      id_start <- merge_ids[[i]][1]
      id_end <- merge_ids[[i]][l]
      lines_list2[[i]] <- c(
        line_coords[id_start, 1:3],
        line_coords[id_end, 4:6]
      )
      df3[i, c("HOLEID", by)] <- df2[id_start, c("HOLEID", by)]
      if (ncol(df3) > 2){
        if (l > 1){ # average
          w <- line_lengths[id_start:id_end] /
            sum(line_lengths[id_start:id_end])
          lweights <- matrix(w, l, ncol(df2) - 2)
          df3[i,seq(3,ncol(df3))] <- colSums(
              lweights * df2[id_start:id_end, seq(3, ncol(df3)), drop=F]
            )
        }else{ # just copy (faster)
          df3[i, seq(3, ncol(df3))] <-
            df2[id_start, seq(3, ncol(df3))]
        }
      }
    }

    # end
    if (by == ".HOLEID") df3 <- df3[, -which(colnames(df3) == by), drop=F]
    return(lines3DDataFrame(lines_list2, df3))
  }
)

#### GetContacts ####
#' @rdname GetContacts
setMethod(
  f = "GetContacts",
  signature = "lines3DDataFrame",
  definition = function(object, by){
    # setup
    x <- MergeSegments(object, by)
    x <- Pointify(x, c(0,1))
    pts <- GetCoords(x,"matrix")
    xdata <- GetData(x)
    # finding duplicate indices
    dup <- which(duplicated(pts))
    dup2 <- dup - 1
    # building new object
    new_points <- GetCoords(x)[dup]
    new_df <- data.frame(xdata[dup2, "HOLEID"], xdata[dup2, by], xdata[dup, by])
    colnames(new_df) <- c("HOLEID", paste0(by, c(".up", ".down")))
    return(points3DDataFrame(new_points, new_df))
  }
)

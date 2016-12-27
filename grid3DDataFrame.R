#### grid3DDataFrame class ####
grid3DDataFrame <- setClass(
        "grid3DDataFrame",
        slots = c(dim = "matrix"),
        contains = "points3DDataFrame",
        validity = function(object) {
                # callNextMethod(object)
        }
)

#### initialization ####
setMethod(
        f = "initialize",
        signature = "grid3DDataFrame",
        definition = function(.Object, gridx, gridy, gridz, fields){
                # dimensions
                nx <- length(gridx)
                ny <- length(gridy)
                nz <- length(gridz)
                .Object@dim <- matrix(c(nx, ny, nz), 1, 3)
                colnames(.Object@dim) <- c("X","Y","Z")
                
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
                p3df <- points3DDataFrame(coords,df)
                .Object@coords <- p3df@coords
                .Object@data <- p3df@data
                .Object@bbox <- p3df@bbox
                # validObject(.Object)
                return(.Object)
        }
)

#### show ####
setMethod(
        f = "show",
        signature = "grid3DDataFrame",
        definition = function(object){
                # setup
                l <- min(10, nrow(object))
                suppressWarnings(p3df <- object[seq(l), ])
                coords <- getCoords(p3df,"data.frame")
                df <- getData(p3df)
                # display
                cat("Object of class ", class(object), "\n", sep = "")
                cat(nrow(object), " coordinates and ",
                    ncol(object), " attributes\n\n", sep = "")
                cat("Number of points:\n")
                d <- object@dim; rownames(d) <- ""; show(d)
                cat("\nBounding Box:\n")
                show(boundingBox(object))
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
        definition = function(x,i,j,drop){
                if(class(i) == "character"){
                        j <- i
                        i <- seq(nrow(x))
                }
                coords_list <- getCoords(x)
                df <- getData(x)
                coords_sub <- coords_list[i]
                df_sub <- df[i,j,drop=FALSE]
                return(points3DDataFrame(coords_sub, df_sub))
                # return(new(class(x), coords_sub, df_sub))
        }
)


#### make3DArray ####
setMethod(
        f = "make3DArray",
        signature = "grid3DDataFrame",
        definition = function(object, value){
                df <- getData(object)
                ar3d <- array(data = df[,value], dim = object@dim)
                coords <- getCoords(object, "matrix")
                return(list(value = ar3d, 
                            x = sort(unique(coords[,1])),
                            y = sort(unique(coords[,2])),
                            z = sort(unique(coords[,3])))
                       )
        }
)

#### reScale ####
# setMethod(
#         f = "reScale",
#         signature = "grid3DDataFrame",
#         definition = function(object, 
#                               old_range = boundingBox(object), 
#                               new_range = matrix(rep(c(0,1),3),2,3)){
#                 points_list <- getCoords(object)
#                 df <- getData(object)
#                 points_list <- lapply(
#                         points_list,
#                         function(x) reScale(x, old_range, new_range)
#                 )
#                 pts <- points3DDataFrame(points_list,df)
#                 object@points <- points_list
#                 object@bbox <- pts@bbox
#                 return(object)
#         }
# )

#### drawSection ####
setMethod(
        f = "drawSection",
        signature = "grid3DDataFrame",
        definition = function(object, by, values, col, 
                              x = NULL, y = NULL, z = NULL){
                # pacakges
                require(rgl)
                require(grDevices)
                require(plotrix)
                # setup
                df <- getData(object)
                if(!any(colnames(df) %in% by)){
                        stop("Invalid attribute")
                }
                ar <- make3DArray(object, by)
                d <- dim(ar$value)
                h <- c(ar$x[2] - ar$x[1], 
                       ar$y[2] - ar$y[1], 
                       ar$z[2] - ar$z[1]) * 0.5
                # pallette
                if(class(as.vector(ar$value)) == "numeric"){ # continuous variable
                        colnum <- col2rgb(col)/255
                        colorsc <- character(length(ar$value))
                        for(i in seq(length(values)-1)){
                                id <- ar$value >= values[i] &
                                        ar$value <= values[i+1]
                                suppressWarnings(
                                        colorsc[id] <- plotrix::color.scale(
                                                ar$value[id],
                                                colnum["red",c(i,i+1)],
                                                colnum["green",c(i,i+1)],
                                                colnum["blue",c(i,i+1)]
                                        )
                                )
                        }
                }else{ # categorical variable
                        for(i in seq_along(values)){
                                ar$value[ar$value == values[i]] <- i
                        }
                        ar$value <- as.numeric(ar$value)
                        colorsc <- col[ar$value]
                        dim(ar$value) <- d
                }
                dim(colorsc) <- d
                # plotting
                if(!is.null(x)){
                        # nearest grid position
                        id <- which.min(abs(x - ar$x)) 
                        # color ordering
                        section <- ar$value[id,,]
                        ord <- sort(section, index.return = T)
                        colorsc <- colorsc[id,,]
                        colorsc <- unique(colorsc[ord$ix])
                        # rescaling section
                        v <- unique(sort(section))
                        for(i in seq_along(v)){
                                section[section == v[i]] <- i
                        }
                        # section
                        show2d({
                                par(mar=c(0,0,0,0))
                                image(section, col = colorsc)
                                grid(d[2], d[3], col = "black", 
                                     lty = "solid")
                        },
                        x = rep(ar$x[id],4),
                        y = c(ar$y[1]-h[2], ar$y[d[2]]+h[2], 
                              ar$y[d[2]]+h[2], ar$y[1]-h[2]),
                        z = c(ar$z[1]-h[3], ar$z[1]-h[3], 
                              ar$z[d[3]]+h[3], ar$z[d[3]]+h[3]),
                        width = 20*d[2], height = 20*d[3]
                        )
                        invisible(NULL)
                }else if(!is.null(y)){
                        # nearest grid position
                        id <- which.min(abs(y - ar$y)) 
                        # color ordering
                        section <- ar$value[,id,]
                        ord <- sort(section, index.return = T)
                        colorsc <- colorsc[,id,]
                        colorsc <- unique(colorsc[ord$ix])
                        # rescaling section
                        v <- unique(sort(section))
                        for(i in seq_along(v)){
                                section[section == v[i]] <- i
                        }
                        # section
                        show2d({
                                par(mar=c(0,0,0,0))
                                image(section, col = colorsc)
                                grid(d[1], d[3], col = "black", 
                                     lty = "solid")
                        },
                        x = c(ar$x[1]-h[1], ar$x[d[1]]+h[1], 
                              ar$x[d[1]]+h[1], ar$x[1]-h[1]),
                        y = rep(ar$y[id],4),
                        z = c(ar$z[1]-h[3], ar$z[1]-h[3], 
                              ar$z[d[3]]+h[3], ar$z[d[3]]+h[3]),
                        width = 20*d[1], height = 20*d[3]
                        )
                        invisible(NULL)
                }else if(!is.null(z)){
                        # nearest grid position
                        id <- which.min(abs(z - ar$z)) 
                        # color ordering
                        section <- ar$value[,,id]
                        ord <- sort(section, index.return = T)
                        colorsc <- colorsc[,,id]
                        colorsc <- unique(colorsc[ord$ix])
                        # rescaling section
                        v <- unique(sort(section))
                        for(i in seq_along(v)){
                                section[section == v[i]] <- i
                        }
                        # section
                        show2d({
                                par(mar=c(0,0,0,0))
                                image(section, col = colorsc)
                                grid(d[1], d[2], col = "black", 
                                     lty = "solid")
                        },
                        x = c(ar$x[1]-h[1], ar$x[d[1]]+h[1], 
                              ar$x[d[1]]+h[1], ar$x[1]-h[1]),
                        y = c(ar$y[1]-h[2], ar$y[1]-h[2], 
                              ar$y[d[2]]+h[2], ar$y[d[2]]+h[2]),
                        z = rep(ar$z[id],4),
                        width = 20*d[1], height = 20*d[2]
                        )
                        invisible(NULL)
                }else{
                        stop("One of x, y, or z must be provided")
                }
        }
)

#### pointify ####
setMethod(
        f = "pointify",
        signature = "grid3DDataFrame",
        definition = function(object){
                points3DDataFrame(getCoords(object), getData(object))
        }
)


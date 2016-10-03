#### spatial3DDataFrame class ####
spatial3DDataFrame <- setClass(
        "spatial3DDataFrame",
        slots = c(coords = "list",
                  data = "data.frame",
                  bbox = "matrix")
)

#### getData ####
setMethod(
        f = "getData",
        signature = "spatial3DDataFrame",
        definition = function(object){
                return(object@data)
        }
)

#### boundingBox ####
setMethod(
        f = "boundingBox",
        signature = "spatial3DDataFrame",
        definition = function(object){
                return(object@bbox)
        }
)

#### nrow, ncol ####
setMethod("nrow", "spatial3DDataFrame",
          function(x){return(nrow(x@data))}
)
setMethod("ncol", "spatial3DDataFrame",
          function(x){return(ncol(x@data))}
)

#### show ####
setMethod(
        f = "show",
        signature = "spatial3DDataFrame",
        definition = function(object){
                # setup
                l <- min(10, nrow(object))
                if(l > 0){
                        sp3df <- object[seq(l), ]
                        coords <- getCoords(sp3df,"data.frame")
                        df <- getData(sp3df)
                        # display
                        cat("Object of class ", class(object), "\n", sep = "")
                        cat(nrow(object), " coordinates and ",
                            ncol(object), " attributes\n\n", sep = "")
                        cat("Bounding Box:\n")
                        show(boundingBox(object))
                        cat("\nCoordinates:\n")
                        show(coords)
                        cat("\nAttributes:\n")
                        show(df)
                } 
                else
                        cat("Empty", class(object))
        }
)

#### [ ####
setMethod(
        f = "[",
        signature = "spatial3DDataFrame",
        definition = function(x,i,j,drop){
                if(class(i) == "character"){
                        j <- i
                        i <- seq(nrow(x))
                }
                coords_list <- getCoords(x)
                df <- getData(x)
                coords_sub <- coords_list[i]
                df_sub <- df[i,j,drop=FALSE]
                # return(points3DDataFrame(points_sub, df_sub))
                return(new(class(x), coords_sub, df_sub))
        }
)

#### [<- ####
setMethod(
        f = "[<-",
        signature = "spatial3DDataFrame",
        definition = function(x, i, j, value, drop){
                if(class(i) == "character"){
                        j <- i
                        i <- seq(nrow(x))
                }
                if(class(value) %in% c("spatial3DDataFrame",
                                       "points3DDataFrame",
                                       "lines3DDataFrame")){
                        value <- getData(value)
                }
                # coords_list <- getCoords(x)
                df <- getData(x)
                df[i,j] <- value
                # return(points3DDataFrame(points_list, df))
                # return(new(class(x), coords_list, df))
                x@data <- df
                return(x)
        }
)

#### [[ ####
setMethod(
        f = "[[",
        signature = "spatial3DDataFrame",
        definition = function(x,i,exact){
                return(unlist(getData(x)[,i]))
        }
)

#### getNormals ####
setMethod(
        f = "getNormals",
        signature = "spatial3DDataFrame",
        definition = function(object, dip = "Dip", strike = "Strike"){
                # setup
                df <- getData(object)
                df[,".diprad"] <- -df[,dip] * pi/180
                df[,".strad"] <- (90 - df[,strike]) * pi/180
                # normal vector
                dipstr <- getPlaneDirections(object, dip, strike)
                vec1 <- getData(dipstr[1:nrow(object), ])
                vec2 <- getData(dipstr[(nrow(object)+1):(2*nrow(object)), ])
                normalvec <- vec1[,c(2,3,1)] * vec2[,c(3,1,2)] -
                        vec1[,c(3,1,2)] * vec2[,c(2,3,1)]
                # normalization
                normalvec <- t(apply(normalvec,1,function(v) v/sqrt(sum(v^2))))
                # result
                colnames(normalvec) <- c("nX", "nY", "nZ")
                normalvec <- as.data.frame(normalvec)
                return(new(class(object), getCoords(object), normalvec))
        }
)

#### getPlaneDirections ####
setMethod(
        f = "getPlaneDirections",
        signature = "spatial3DDataFrame",
        definition = function(object, dip = "Dip", strike = "Strike"){
                # setup
                df <- getData(object)
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
                coords <- getCoords(object)
                coords <- c(coords, coords)
                return(new(class(object), coords, vecs))
        }
)

#### getLineDirections ####
setMethod(
        f = "getLineDirections",
        signature = "spatial3DDataFrame",
        definition = function(object, dip = "Dip", azimuth = "Azimuth"){
                # setup
                df <- getData(object)
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
                coords <- getCoords(object)
                return(new(class(object), coords, dipvec))
        }
)
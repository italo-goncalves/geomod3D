#### line3D class ####
line3D <- setClass("line3D",
                   slots = c(coords = "matrix",
                             length = "numeric"),
                   prototype = list(coords = matrix(c(0,0,0,1,1,1),
                                                    nrow = 2,
                                                    ncol = 3,
                                                    byrow = TRUE),
                                    length = sqrt(3)),
                   validity = function(object){
                           if(!all(dim(object@coords) == c(2,3))){
                                   stop(paste("Vectors with length",
                                              "of 2 or 3 must be provided"))
                           }
                           if(!is.numeric(object@coords)){
                                   stop("coords must be numeric")
                           }
                           return(TRUE)
                   }
)

#### initialization ####
setMethod(
        f = "initialize",
        signature = "line3D",
        definition = function(.Object, pointA, pointB){
                if(class(pointA) == "point3D"){
                        pointA <- getCoords(pointA)
                }
                if(class(pointB) == "point3D"){
                        pointB <- getCoords(pointB)
                }
                if(class(pointA) == "data.frame"){
                        pointA <- as.matrix(pointA)
                }
                if(class(pointB) == "data.frame"){
                        pointB <- as.matrix(pointB)
                }
                if(length(pointA) == 2){
                        pointA <- append(pointA, 0)
                }
                if(length(pointB) == 2){
                        pointB <- append(pointB, 0)
                }
                .Object@coords <- matrix(c(pointA,pointB),
                                         nrow = 2,
                                         ncol = 3,
                                         byrow = TRUE)
                .Object@length <- sqrt(sum((.Object@coords[1,] - 
                                                    .Object@coords[2,])^2))
                colnames(.Object@coords) <- c("X","Y","Z")
                validObject(.Object)
                return(.Object)
        }
)

#### getCoords ####
setMethod(
        f = "getCoords",
        signature = "line3D",
        definition = function(object){
                return(object@coords)
        }
)

#### getLength ####
setMethod(
        f = "getLength",
        signature = "line3D",
        definition = function(object){
                return(object@length)
        }
)

#### discretize ####
# setMethod(
#         f = "discretize",
#         signature = "line3D",
#         definition = function(object, interval, drop.stubs = TRUE){
#                 # fazer
#         }
# )

#### pointify ####
setMethod(
        f = "pointify",
        signature = "line3D",
        definition = function(object, locations = c(0.05,0.5,0.95)){
                # validation
                if(min(locations) < 0 | max(locations) > 1) {
                        stop(paste("argument 'locations' must be a vector",
                                   "with values between 0 and 1"))
                }
                # conversion to point3D
                points_list <- vector("list", length(locations))
                for(i in seq_along(locations)){
                        coords <- getCoords(object)
                        points_list[[i]] <- point3D(
                                coords[1,] + locations[i] * 
                                        (coords[2,]-coords[1,])
                        )
                        validObject(points_list[[i]])
                }
                return(points_list)
        }
)

#### reScale ####
setMethod(
        f = "reScale",
        signature = "line3D",
        definition = function(object, 
                              old_range = getCoords(object), 
                              new_range = matrix(rep(c(0,1),3),2,3)){
                coords <- getCoords(object)
                coords <- rescale_coords(coords, old_range, new_range)
                return(line3D(coords[1,],coords[2,]))
        }
)
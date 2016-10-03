#### point3D class ####
point3D <- setClass("point3D",
                    slots = c(coords = "matrix"),
                    prototype = list(coords = matrix(numeric(3),
                                                     nrow = 1,
                                                     ncol = 3)),
                    validity = function(object){
                            if(length(object@coords) != 3){
                                    stop(paste("A vector with length",
                                               "of 2 or 3 must be provided"))
                            }
                            if(!is.numeric(object@coords)){
                                    stop("coords must be numeric")
                            }
                            return(TRUE)
                    }
)

#### initalization ####
setMethod(
        f = "initialize",
        signature = "point3D",
        definition = function(.Object, coords){
                if(class(coords) == "data.frame"){
                        coords <- as.matrix(coords)
                }
                if(length(coords) == 2){
                        coords <- append(coords, 0)
                }
                .Object@coords <- matrix(coords,
                                         nrow = 1,
                                         ncol = 3)
                colnames(.Object@coords) <- c("X","Y","Z")
                # validObject(.Object)
                return(.Object)
        }
)

#### getCoords ####
setMethod(
        f = "getCoords",
        signature = "point3D",
        definition = function(object){
                return(object@coords)
        }
)

#### reScale ####
setMethod(
        f = "reScale",
        signature = "point3D",
        definition = function(object, old_range, new_range){
                coords <- getCoords(object)
                coords <- rescale_coords(coords, old_range, new_range)
                object@coords <- matrix(coords,1,3)
                return(object)
        }
)
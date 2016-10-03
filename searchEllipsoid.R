#### searchEllipsoid class ####
searchEllipsoid <- setClass(
        "searchEllipsoid",
        slots = c(
                maxrange = "numeric",
                midrange = "numeric",
                minrange = "numeric",
                azimuth = "numeric",
                dip = "numeric",
                rake = "numeric",
                sectors = "numeric",
                points.min = "numeric",
                points.max = "numeric"
                ),
        validity = function(object) {
                at <- list(object@maxrange, object@midrange,
                           object@minrange, object@azimuth, object@dip,
                           object@rake, object@sectors, object@points.min,
                           object@points.max)
                if(any(lapply(at, class) != "numeric"))
                        stop("All ellipsoid parameters must be numeric
                             values")
                return(TRUE)
        }
)

#### initialization ####
setMethod(
        f = "initialize",
        signature = "searchEllipsoid",
        definition = function(.Object, maxrange, midrange = maxrange,
                              minrange = midrange, azimuth = 0, dip = 0,
                              rake = 0, sectors = 4, points.min = 1,
                              points.max = 5){
                .Object@maxrange = maxrange
                .Object@midrange = midrange
                .Object@minrange = minrange
                .Object@azimuth = azimuth
                .Object@dip = dip
                .Object@rake = rake
                .Object@sectors = sectors
                .Object@points.min = points.min
                .Object@points.max = points.max
                validObject(.Object)
                return(.Object)
        }
)

#### show ####
setMethod(
        f = "show",
        signature = "searchEllipsoid",
        definition = function(object){
                cat("Object of class ", class(object), "\n\n", sep = "")
                cat("Maximum range =", object@maxrange,"\n")
                cat("Intermediate range =", object@midrange,"\n")
                cat("Minimum range =", object@minrange,"\n")
                cat("Orientation: azimuth =", object@azimuth,
                    "dip =", object@dip, "rake =", object@rake, "\n")
                cat("Number of setors =", object@sectors,"\n")
                cat("Minimum number of points per sector =", 
                    object@points.min,"\n")
                cat("Maximum number of points per sector =", 
                    object@points.max,"\n")
        }
)
############### generic functions #############################

#### getters ####
setGeneric("getCoords",
           function(object, ...){standardGeneric("getCoords")}
)

setGeneric("getLength",
           function(object, ...){standardGeneric("getLength")}
)

# setGeneric("getLines",
#            function(object, ...){standardGeneric("getLines")}
# )
# 
# setGeneric("getPoints",
#            function(object, ...){standardGeneric("getPoints")}
# )

setGeneric("getData",
           function(object, ...){standardGeneric("getData")}
)

setGeneric("boundingBox",
           function(object, ...){standardGeneric("boundingBox")}
)

#### specialized functions ####
# setGeneric("discretize",
#            function(object, ...){standardGeneric("discretize")}
# )

setGeneric("pointify",
           function(object, ...){standardGeneric("pointify")}
)

# setGeneric("reScale",
#            function(object, ...){standardGeneric("reScale")}
# )

setGeneric("make3DArray",
           function(object, ...){standardGeneric("make3DArray")}
)

setGeneric("mergeSegments",
           function(object, ...){standardGeneric("mergeSegments")}
)

setGeneric("getContacts",
           function(object, ...){standardGeneric("getContacts")}
)

setGeneric("bindPoints",
           function(x, y, ...){standardGeneric("bindPoints")}
)

#### structural geology ####
setGeneric("getNormals",
           function(object, ...){standardGeneric("getNormals")}
)

setGeneric("getPlaneDirections",
           function(object, ...){standardGeneric("getPlaneDirections")}
)

setGeneric("getLineDirections",
           function(object, ...){standardGeneric("getLineDirections")}
)

setGeneric("contactsBestFitPlane",
           function(object, ...){standardGeneric("contactsBestFitPlane")}
)

#### geostatistics ####
setGeneric("covariance_matrix",
           function(x, y, ...){standardGeneric("covariance_matrix")}
)

setGeneric("covariance_matrix_d1",
           function(x, tangents, ...){standardGeneric("covariance_matrix_d1")}
)

setGeneric("covariance_matrix_d2",
           function(tangents, ...){standardGeneric("covariance_matrix_d2")}
)

setGeneric("trend_matrix",
           function(x, trend, ...){standardGeneric("trend_matrix")}
)

setGeneric("trend_matrix_d1",
           function(x, trend, ...){standardGeneric("trend_matrix_d1")}
)

# setGeneric("krig3D",
#            function(x, y, ...){standardGeneric("krig3D")}
# )
# 
# setGeneric("krig3D_search",
#            function(x, y, ...){standardGeneric("krig3D_search")}
# )
# 
# setGeneric("krig3D_cluster",
#            function(x, y, ...){standardGeneric("krig3D_cluster")}
# )

# setGeneric("smoothingPenalty",
#            function(object, ...){standardGeneric("smoothingPenalty")}
# )

# setGeneric("classify",
#            function(x, y, ...){standardGeneric("classify")}
# )

# setGeneric("classify_ps",
#            function(x, y, ...){standardGeneric("classify_ps")}
# )


#### visualization ####
setGeneric("drawDrillholes",
           function(object, ...){standardGeneric("drawDrillholes")}
)

setGeneric("drawSection",
           function(object, ...){standardGeneric("drawSection")}
)

setGeneric("drawPoints",
           function(object, ...){standardGeneric("drawPoints")}
)

setGeneric("drawTangentPlanes",
           function(object, ...){standardGeneric("drawTangentPlanes")}
)

setGeneric("drawTangentLines",
           function(object, ...){standardGeneric("drawTangentLines")}
)

setGeneric("drawHoleID",
           function(object, ...){standardGeneric("drawHoleID")}
)

#### GP objects ####
setGeneric("fit",
           function(object, ...){standardGeneric("fit")}
)
setGeneric("fit_gating",
           function(object, ...){standardGeneric("fit_gating")}
)
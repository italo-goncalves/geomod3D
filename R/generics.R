############### generic functions #############################
#' @include geomod3D.R
#' @include utils.R
NULL

#### getters ####

#' \code{spatial3DDataFrame} objects' coordinates
#'
#' Returns the spatial coordinates of a \code{spatial3DDataFrame} object.
#'
#' @param object The object to get coordinates from.
#' @param as The return format. One of \code{"list"}, \code{"matrix"}, or
#' \code{"data.frame"}.
#'
#' @return An object according to the \code{as} parameter containing the
#' object's coordinates.
#'
#' @export
setGeneric("GetCoords",
           function(object, ...){standardGeneric("GetCoords")}
)

#' Drillhole segments' lengths
#'
#' Returns the length of each segment of a \code{lines3DDataFrame} object.
#'
#' @param object The object to get the data from.
#' @return A numeric vector with the segments' lengths.
#' @export
setGeneric("GetLength",
           function(object, ...){standardGeneric("GetLength")}
)

#' \code{spatial3DDataFrame} objects' data
#'
#' Returns the attributes of a \code{spatial3DDataFrame} object, represented
#' as a \code{data.frame}.
#'
#' @param object The object to get the data from.
#' @return A data frame with the object's data.
#' @export
setGeneric("GetData",
           function(object, ...){standardGeneric("GetData")}
)

#' Bounding box
#'
#' Returns a \code{spatial3DDataFrame} object's, bounding box, the smallest
#' rectangular parallelepiped that contains all the data.
#'
#' @param object The object from which to obtain the bounding box.
#'
#' @return A 2x3 matrix containing the lower and upper bounds of the object's
#' coordinates.
#' @export
setGeneric("BoundingBox",
           function(object, ...){standardGeneric("BoundingBox")}
)

#### specialized functions ####
#' Conversion to point cloud
#'
#' Converts any \code{spatial3DDataFrame} object to a \code{points3DDataFrame}.
#' Specific options depend of the object being converted.
#' @param object The object being converted.
#' @param locations Relative locations along each line segment from which to
#' extract a point (between 0 and 1, inclusive).
#' @param distance Should the distance from the closest vertex be included?
#'
#' @details If \code{distance = T}, the distance from the point to the closest
#' vertex in the original line segment is included in the returned object. The
#' values are stored in the \code{.dist} column.
#'
#' In the current version, this method is only useful for
#' \code{lines3DDataFrame} objects. Future versions will add functionality for
#' block models and other kinds of objects.
#'
#' @return A \code{points3DDataFrame} object.
setGeneric("Pointify",
           function(object, ...){standardGeneric("Pointify")}
)

#' MergeSegments
#'
#' Merging of redundant line segments to reduce the number of rows in a
#' \code{lines3DDataFrame} object.
#'
#' @param by Column names to be used on merging.
#' @param keep Which columns to keep in the merged object.
#'
#' @details This method merges contiguous line segments that are parallel and
#' for which all variables described in the \code{by} parameter are equal.
#' Variables described in the \code{keep} parameter will be averaged if numeric.
#' If any variable described in \code{keep} is a character of factor, new
#' columns will be created in the merged object that contain the relative
#' proportion of each unique element or factor in each merged segment.
#'
#' @return A \code{lines3DDataFrame} object.
setGeneric("MergeSegments",
           function(object, ...){standardGeneric("MergeSegments")}
)

#' GetContacts
#'
#' Finds the points of contact between different geological units or other
#' categorical variable.
#'
#' @param by A column name indicating which attribute to use to find contacts.
#'
#' @details This method extracts the points that lie in the border between any
#' two contiguous line segments with different values of the \code{by} column.
#' It is assumed that \code{by} represents a categorical variable.
#'
#' @return A \code{points3DDataFrame} object with three attributes:
#' \describe{
#'   \item{HOLEID}{The identification of the hole from which the point was
#'   extracted.}
#'   \item{.up}{The value of the \code{by} variable in the segment above the
#'   point}
#'   \item{.down}{The value of the \code{by} variable in the segment below the
#'   point}
#' }
setGeneric("GetContacts",
           function(object, ...){standardGeneric("GetContacts")}
)

#' Bind
#'
#' Binds two 3D spatial objects together.
#'
#' @param x,y \code{spatial3DDataFrame} objects.
#'
#' @details The two objects' coordinates are binded as in the \code{rbind}
#' function. If the object's column names are not the same, the resulting
#' object will have all columns from both, with the appropriate entries filled
#' with \code{NA}. Note that if the parent objects have coordinates in common,
#' the resulting object will have duplicated coordinates.
#'
#' @return The resulting object's class will depend of the classes of the
#' parent objects.
setGeneric("Bind",
           function(x, y, ...){standardGeneric("Bind")}
)

#' Make3DArray
#'
#' Converts a regularly spaced 3D object to a 3D array.
#'
#' @param value The column name to use. It can represent either a continuous or
#' a categorical variable.
#'
#' @return A \code{list} containing the following:
#' \describe{
#'   \item{value}{The 3D array itself.}
#'   \item{x, y, z}{The coordinates in each direction.}
#' }
setGeneric("Make3DArray",
           function(object, ...){standardGeneric("Make3DArray")}
)

#### structural geology ####
#' Normals to structural planes
#'
#' Calculates the vector normal to a structural plane.
#'
#' @param object A \code{spatial3DDataFrame} object.
#' @param dip,strike Names of columns containing dip and strike information.
#'
#' @details Dip and strike values are assumed to be in degrees.
#'
#' @return An object of the same class as \code{object} containing the
#' coordinates of the normals to the specified planes (nX, nY, nZ). The
#' coordinates are normalized to unit length.
setGeneric("GetNormals",
           function(object, ...){standardGeneric("GetNormals")}
)

#' Structural planes' directions
#'
#' Calculates structural planes' dip and strike vectors.
#'
#' @param object A \code{spatial3DDataFrame} object.
#' @param dip,strike Names of columns containing dip and strike information.
#'
#' @details Dip and strike values are assumed to be in degrees.
#'
#' @return An object of the same class as \code{object} (with twice the number
#' of rows) containing the coordinates of the vectors that define a structural
#' plane (one parallel to dip and the other parallel to strike). The
#' coordinates are normalized to unit length.
setGeneric("GetPlaneDirections",
           function(object, ...){standardGeneric("GetPlaneDirections")}
)

#' Structural lines' directions
#'
#' Calculates the coordinates of vectors corresponding to structural directions.
#'
#' @param object A \code{spatial3DDataFrame} object.
#' @param dip,azimuth Names of columns containing dip and azimuth information.
#'
#' @details Dip and azimuth values are assumed to be in degrees.
#'
#' @return An object of the same class as \code{object} containing the
#' coordinates of the vectors that define a structural
#' line. The coordinates are normalized to unit length.
setGeneric("GetLineDirections",
           function(object, ...){standardGeneric("GetLineDirections")}
)


#### geostatistics ####
#' Covariance matrix of spatial data and/or its derivatives.
#'
#' Calculates the covariance matrix of a spatially distributed value and/or its
#' derivative.
#'
#' @param x,y 3D spatial objects.
#' @param tangents A \code{points3DDataFrame} containing directional data, most
#' likely generated with the \code{GetPlaneDirections()} method.
#' @param model A \code{covarianceStructure3D} object representing the spatial
#' continuity of the variable of interest, or a \code{list} containing multiple
#' such objects.
#' @param covariance Should the matrix be calculated as a covariance matrix
#' (default) or a variance (gamma) matrix?
#' @param parts The number of parts in which to "break" line segments.
#'
#' @details If \code{covariance = F} the resulting matrix is given in variogram
#' form, with zeroes in the main diagonal and values that increase with the
#' distance between data locations.
#'
#' If the spatial objects are not made of punctual data, the covariance
#' is calcuated with a "brute force" approach, by discretizing each data
#' element in points and calculating an average covariance.
#'
#' The \code{CovarianceMatrixD1()} method calculates the covariances between
#' the variable of interest located in \code{x} and its derivatives located in
#' \code{tangents}. The \code{CovarianceMatrixD2()} calculates the covariances
#' between the derivatives themselves.
#'
#' @return A covariance matrix with dimensions depending of the method called:
#' \describe{
#'   \item{\code{CovarianceMatrix()}}{ \code{nrow(x)} rows and \code{nrow(y)}
#'   columns}
#'   \item{\code{CovarianceMatrixD1()}}{ \code{nrow(x)} rows and
#'   \code{nrow(tangents)} columns}
#'   \item{\code{CovarianceMatrixD2()}}{ \code{nrow(tangents)} rows and
#'   \code{nrow(tangents)} columns}
#' }
#'
#' @seealso \code{\link{covarianceStructure3D-class}},
#' \code{\link{GetPlaneDirections}}
setGeneric("CovarianceMatrix",
           function(x, y, ...){standardGeneric("CovarianceMatrix")}
)

#' @rdname CovarianceMatrix
setGeneric("CovarianceMatrixD1",
           function(x, tangents, ...){standardGeneric("CovarianceMatrixD1")}
)

#' @rdname CovarianceMatrix
setGeneric("CovarianceMatrixD2",
           function(tangents, ...){standardGeneric("CovarianceMatrixD2")}
)

#' Trend matrix
#'
#' Calculates a matrix containing the components of a trend function at the
#' specified data locations.
#'
#' @param x A 3D spatial object.
#' @param trend A character string convertable to a formula.
#'
#' @details The \code{TrendMatrixD1()} method calculates the derivative of the
#' trend function at the given locations. This is used in kriging systems
#' involving derivative data.
#'
#' @return A matrix with \code{nrow(x)} rows and number of columns according to
#' the specified trend formula.
#'
#' @seealso \code{\link{CovarianceMatrix}}
setGeneric("TrendMatrix",
           function(x, trend, ...){standardGeneric("TrendMatrix")}
)

#' @rdname TrendMatrix
setGeneric("TrendMatrixD1",
           function(x, trend, ...){standardGeneric("TrendMatrixD1")}
)

#### visualization ####
#' Visualization of drillhole data
#'
#' Adds drillhole information to the current \code{rgl} window as colored
#' cylinders.
#'
#' @param by Name of the column to show.
#' @param values,col If \code{by} is a categorical variable, the values to show
#' and the corresponding color. If \code{by} is a continuous variable, a numeric
#' vector with the breakpoints and the corresponding colors.
#' @param size The diameter of the cylinders. Either a single value or a vector
#' matching the object's number of rows.
#' @param col.default Color to be used for \code{NA} values and values outside
#' the range provided.
#'
#' @details For categorical variables, each entry in \code{values} is matched to
#' the corresponding entry in \code{col}. For continuous variables,
#' \code{values} is a numeric vector with breakpoints, each of which is matched
#' to the entries in \code{col}. The actual values are interpolated to generate
#' a unique color to each value. \code{col} must be a character vector
#' containing colors in hexadecimal format, such as \code{"#804DB3"}, or valid
#' color names.
setGeneric("DrawDrillholes",
           function(object, ...){standardGeneric("DrawDrillholes")}
)

#' Visualization of drillhole data
#'
#' Draws the hole identifications in the current \code{rgl} window.
#'
#' @param cex The character expansion factor for the labels.
setGeneric("DrawHoleID",
           function(object, ...){standardGeneric("DrawHoleID")}
)

#' Visualization of point data
#'
#' Adds point information to the current \code{rgl} window in the form of
#' colored spheres or points.
#'
#' @param by Name of the column to show.
#' @param values,col If \code{by} is a categorical variable, the values to show
#' and the corresponding color. If \code{by} is a continuous variable, a numeric
#' vector with the breakpoints and the corresponding colors.
#' @param size The diameter of the spheres. Either a single value or a vector
#' matching the object's number of rows.
#' @param alpha Opacity factor.
#' @param col.default Color to be used for \code{NA} values and values outside
#' the range provided.
#' @param as How to plot the points. For \code{"points"}, only the color
#' information is used but the visualization becomes much lighter in terms of
#' machine resources.
#'
#' @details For categorical variables, each entry in \code{values} is matched to
#' the corresponding entry in \code{col}. For continuous variables,
#' \code{values} is a numeric vector with breakpoints, each of which is matched
#' to the entries in \code{col}. The actual values are interpolated to generate
#' a unique color to each value. \code{col} must be a character vector
#' containing colors in hexadecimal format, such as \code{"#804DB3"}, or valid
#' color names.
setGeneric("DrawPoints",
           function(object, ...){standardGeneric("DrawPoints")}
)

#' Visualization of structural data
#'
#' Draws discs representing the orientation of structural planes.
#'
#' @param size The diameter of the discs.
#' @param dip,strike Names of the columns that contain the structural
#' information.
#' @param col The color of the discs. Either a single value or a vector
#' matching the object's number of rows.
#'
#' @details \code{col} must contain colors in hexadecimal format, such as
#' \code{"#804DB3"}, or valid color names.
#'
#' It is assumed that the values in the \code{dip} and \code{strike} columns
#' are in degrees.
setGeneric("DrawTangentPlanes",
           function(object, ...){standardGeneric("DrawTangentPlanes")}
)

#' Visualization of structural data
#'
#' Draws thin cylinders representing the orientation of structural lines.
#'
#' @param size The length of the cylinders.
#' @param dX,dY,dZ Names of the columns that contain the structural lines'
#' directions.
#' @param col The color of the discs. Either a single value or a vector
#' matching the object's number of rows.
#'
#' @details \code{col} must contain colors in hexadecimal format, such as
#' \code{"#804DB3"}, or valid color names.
#'
#' It is assumed that \code{dX}, \code{dY}, and \code{dZ} are coordinates of
#' a unit vector.
setGeneric("DrawTangentLines",
           function(object, ...){standardGeneric("DrawTangentLines")}
)

#' DrawSection
#'
#' Draws a 2D section of a regularly spaced 3D object in the current \code{rgl}
#' window.
#'
#' @param by Name of the column to show.
#' @param values,col If \code{by} is a categorical variable, the values to show
#' and the corresponding color. If \code{by} is a continuous variable, a numeric
#' vector with the breakpoints and the corresponding colors.
#' @param col.default Color to be used for \code{NA} values and values outside
#' the range provided.
#' @param x,y,z The coordinate to be kept fixed in order to draw the section (
#' will be matched to the closest corresponding coordinate in \code{object}).
#' Only one of these must be provided.
#'
#' @details The section is drawn by keeping one coordinate fixed. For example,
#' if \code{x=5} is given, the section will be drawn in the yz plane.
setGeneric("DrawSection",
           function(object, ...){standardGeneric("DrawSection")}
)

#### GP objects ####
#' Predict
#'
#' Spatial interpolation at new locations.
#'
#' @param target The \code{spatial3DDataFrame} object to receive the prediction.
#' @param to The name of the column in which to write the prediction. Will be
#' overwritten if it exists or created if not.
#' @param output.var Should the predictive variance be computed?
#'
#' @return A 3D spatial object of the same class as \code{target}
#' containing the predictions.
setGeneric("Predict",
           function(object, ...){standardGeneric("Predict")}
)

#' Fit
#'
#' Training of a Gaussian Process model with a genetic algorithm.
#'
#' @param midrange,minrange Optimize on the covariance model's anisotropy?
#' @param azimuth,dip,rake Optimize on the covariance model's orientation?
#' @param nugget Optimize on nugget?
#' @param power Optimize on the power parameter (not relevant for all
#' covariances)?
#'
#' @details This method uses a genetic algorithm to optimize the contribution
#' (or amplitude) and range of each covariance structure contained in the
#' \code{model} slot of \code{object}, as well as the parameters above, if
#' allowed. Optimization is done with respect to the object's log-likelihood.
#'
#' By default the genetic algorithm will run for a maximum of 100 iterations or
#' 20 iterations without improvement of the log-likelihood. If desired, this
#' method can be called multiple times to continue the optimization starting
#' from previous results.
#'
#' In order to obtain reproducible results, set the seed before calling this
#' method.
#'
#' @return A \code{GP} object similar to \code{object}, with optimized
#' covariance parameters.
setGeneric("Fit",
           function(object, ...){standardGeneric("Fit")}
)

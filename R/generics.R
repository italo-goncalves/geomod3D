############### generic functions #############################
#' @include geomod3D.R
#' @include utils.R
NULL

#### base ####
setGeneric("ncol")
setGeneric("nrow")
setGeneric("summary")
setGeneric("as.data.frame")
setGeneric("str")

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
#' @return A \code{directions3DDataFrame} containing the
#' coordinates of the normals to the specified planes. The
#' coordinates are normalized to unit length.
#' @seealso \code{\link{GetPlaneDirections}}, \code{\link{GetLineDirections}}
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
#' @return A \code{directions3DDataFrame} (with twice the number of rows as the
#' input) containing the coordinates of the vectors that define a structural
#' plane (one parallel to dip and the other parallel to strike). The
#' coordinates are normalized to unit length.
#' @seealso \code{\link{GetNormals}}, \code{\link{GetLineDirections}}
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
#' @return A \code{directions3DDataFrame} containing the
#' coordinates of the vectors that define a structural
#' line. The coordinates are normalized to unit length.
#' @seealso \code{\link{GetPlaneDirections}}, \code{\link{GetNormals}}
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
#' @param model A \code{covarianceStructure3D} object representing the spatial
#' continuity of the variable of interest, or a \code{list} containing multiple
#' such objects.
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
           function(x, y, model, ...){standardGeneric("CovarianceMatrix")}
)

setGeneric("NuggetMatrix",
           function(x, model, ...){standardGeneric("NuggetMatrix")}
)

#' Trend matrix
#'
#' Calculates a matrix containing the components of a trend function at the
#' specified data locations.
#'
#' @param x A 3D spatial object.
#' @param trend A character string convertable to a formula.
#'
#' @details A call with a \code{directions3DDataFrame} as argument returns
#' the derivative of the trend function at the given locations. This is used in
#' situations involving derivative data.
#'
#' @return A matrix with \code{nrow(x)} rows and number of columns according to
#' the specified trend formula.
#'
#' @seealso \code{\link{CovarianceMatrix}}
setGeneric("TrendMatrix",
           function(x, trend, ...){standardGeneric("TrendMatrix")}
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
#' Draws thin cylinders representing the orientation of structural lines or
#' directional data in general.
#'
#' @param size The length of the cylinders.
#' @param col The color of the discs. Either a single value or a vector
#' matching the object's number of rows.
#' @param as How to plot the data. \code{"direction"} is useful for
#' structural data, while \code{"arrow"} can be used for for velocities and
#' derivative data in a more strict sense.
#'
#' @details \code{col} must contain colors in hexadecimal format, such as
#' \code{"#804DB3"}, or valid color names.
setGeneric("DrawDirections",
           function(object, ...){standardGeneric("DrawDirections")}
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
#' Predictions at new spatial locations.
#'
#' @param target The \code{spatial3DDataFrame} object to receive the prediction.
#' @param to The name of the column in which to write the prediction. Will be
#' overwritten if it exists or created if not.
#' @param output.var Should the predictive variance be computed?
#' @param output.ind Return indicators for boundary drawing?
#' @param output.prob Return class probabilities?
#' @param use.unknown Include the unknown class in output?
#' @param Nsamp Number of samples used to estimate class probabilities.
#'
#' @return A 3D spatial object of the same class as \code{target}
#' containing the predictions.
#'
#' @details \code{GP} and \code{SPGP} objects return the predicted mean and
#' variance for each location in \code{target}. The sparse GP returns two
#' variances: \code{var_full} represents the total prediction uncertainty while
#' \code{var_cor} is the amplitude of variation of the underlying latent
#' function. The proportion between the two varies according to the distance
#' from the pseudo_inputs.
#'
#' The \code{GP_geomod} object will calculate an indicator and its variance
#' for each class at each location, which jointly form a multivariate normal
#' distribution of the true indicators (or log-transformed compositional
#' coordinates). The returned probabilities are actually the proportion of this
#' probability mass over the region in which each indicator is dominant. This
#' quantity is approximated by drawing \code{Nsamp} samples from the
#' distribution and computing the number of times each indicator is higher
#' than the others. The \code{use.unknown = F} option simply drops the
#' probability of the unknown class and re-normalizes the rest in order for
#' them to add to 1.
setGeneric("Predict",
           function(object, ...){standardGeneric("Predict")}
)

#' Fit
#'
#' Training of a Gaussian Process model with a genetic algorithm.
#'
#' @param contribution Optimize on the covariance model's amplitude?
#' @param maxrange Optimize on the covariance model's range?
#' @param midrange,minrange Optimize on the covariance model's anisotropy?
#' @param azimuth,dip,rake Optimize on the covariance model's orientation?
#' @param nugget Optimize on nugget?
#' @param power Optimize on the power parameter (not relevant for all
#' covariances)?
#' @param pseudo_inputs Optimize on the pseudo-inputs' locations?
#' @param pseudo_tangents Optimize on the pseudo-inputs' locations for
#' derivative data?
#' @param metric Which metric to optimize?
#' @param ... Arguments passed on to \code{\link{ga}}, such as \code{maxiter},
#' \code{popSize}, etc.
#'
#' @details This method uses a genetic algorithm to optimize the contribution
#' (or amplitude) and range of each covariance structure contained in the
#' \code{model} slot of \code{object}, as well as the parameters above, if
#' allowed.
#'
#' Optimization is done with respect to the specified metric. The available
#' metrics are the log-likelihood (default), normalized root mean square error
#' (NRMSE), and penalized log predictive density (PLPD). The latter two are
#' determined by leave-one-out cross validation, which is slower than the
#' log-likelihood.
#'
#' The positions of the pseudo-inputs may be constrained to match a subset of
#' the data (\code{pseudo_inputs = "subset"}) or be free to lie anywhere
#' inside the data`s bounding box (\code{pseudo_inputs = "free"}). The
#' pseudo-tangents, however, must be a subset of the tangent data to avoid
#' overfitting, as the tangents have a greater degree of freedom (position and
#' direction).
#'
#' See the documentation in \code{\link{ga}} to control the optimization
#' process. Standard GP uses continuous optimization to fit the parameters,
#' using the current ones as a starting point. The sparse GP uses discrete
#' optimization (with binary encoding) to fit the parameters and select the
#' best positions for the pseudo-inputs. In both cases, it is
#' recommnded to set the \code{popSize} around 20 and
#' \code{pmutation} between 0.3 and 0.5. Convergence status can be visualized by
#' setting \code{monitor = T}.
#'
#' The variational SPGP model may pose difficulties for training. It may help
#' to train a FIC model (or a standard GP, using a subset of the data)
#' to obtain the best covariance
#' parameters and nugget and use them to build a variational SPGP object.
#'
#' In order to obtain reproducible results, set the seed before calling this
#' method.
#'
#' @return A \code{GP} object similar to \code{object}, with optimized
#' covariance parameters.
#'
#' @seealso \code{\link{covarianceStructure3D-class}},
#' \code{\link{GP-init}}, \code{\link{SPGP-init}}
#'
#' @references
#' Bauer, M.S., van der Wilk, M., Rasmussen, C.E., 2016. Understanding
#' Probabilistic Sparse Gaussian Process Approximations. Adv. Neural Inf.
#' Process. Syst. 29.
setGeneric("Fit",
           function(object, ...){standardGeneric("Fit")}
)

#' Geostatistical simulation
#'
#' Generates a number of equally possible realizations of the modeled random
#' field.
#'
#' @param object A \code{GP} object.
#' @param target The \code{spatial3DDataFrame} object to receive the prediction.
#' @param Nsim The desired number of simulations.
#' @param to The name of the column in which to write the prediction. Will be
#' overwritten if it exists or created if not.
#' @param discount.noise Whether to sample from the GP latent function instead
#' of the noisy values.
#' @param smooth Whether to correct the output from the sparse GP.
#' @param verbose Display status while running?
#'
#' @details Standard GP objects use the Cholesky method, which limits the number
#' of data to be used and the number of test points in which to simulate. SPGP
#' objects use a sequential method based on rank-one updates to the matrices
#' involved, which only limits the number of pseudo-inputs that can be used.
#'
#' Due to the nature of the sparse approximation of the covariance matrices,
#' the output may appear noisy even if \code{discount.noise = T}. Set
#' \code{smooth = T} to apply a correction to the generated samples.
#'
#' @return A 3D spatial object of the same class as \code{target}
#' containing the simulations.
#'
setGeneric("Simulate",
           function(object, ...){standardGeneric("Simulate")}
)

#' Cross-validation
#'
#' Performs cross-validation and calculates a number of performance metrics
#'
#' @param object A \code{SPGP} object.
#'
#' @details The method performs leave-one-out cross-validation through
#' rank-one downdates to a covariance matrix.
#'
#' @return A \code{list} with the following elements:
#' \describe{
#'   \item{mean}{The predictive mean for each data point.}
#'   \item{var}{The predictive variance.}
#'   \item{RMSE}{Root Mean Square Error of prediction, taking only the mean in
#'   consideration}
#'   \item{NRMSE}{Normalized RMSE: the squared difference is divided by the
#'   predictive variance.}
#'   \item{LPD}{Log Predictive Density.}
#'   \item{PLPD}{Penalized LPD: includes a penalty based on how well the full
#'   covariance is approximated at each data point.}
#' }
setGeneric("Xval",
           function(object, ...){standardGeneric("Xval")}
)

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GP.R
\name{GP-init}
\alias{GP-init}
\alias{GP}
\title{Gaussian Process}
\usage{
GP(data, model, value, mean = NULL, trend = NULL,
  force.interp = numeric(), reg.v = 1e-09, tangents = NULL,
  reg.t = 1e-12, nugget.t = 0)
}
\arguments{
\item{data}{A \code{spatial3DDataFrame} object containing the data one
wishes to model.}

\item{model}{The covariance model. A \code{covarianceModel3D} object.}

\item{value}{The column name of the variable to be modeled. It is assumed
the column does not contain missing values.}

\item{mean}{The global mean. Irrelevant if a trend is provided.}

\item{trend}{The model's trend component. A formula in character format.}

\item{force.interp}{Indices of points that must be interpolated exactly.}

\item{reg.v}{Regularization to improve stability. A single value or a vector
with length matching the number of data points.}

\item{tangents}{A \code{directions3DDataFrame} object containing structural
geology data. Most likely generated with the \code{GetLineDirections()}
method.}

\item{reg.t}{Regularization for structural data. A single value or a vector
with length matching the number of structural data.}

\item{weights}{The importance of each data point in the model (a vector with
values between 0 and 1)}
}
\description{
Implementation of the Gaussian Process model for 3D spatial interpolation.
}
\details{
This method builds a \code{GP} object with all the information
needed to make preditions at new data points.

\code{trend} must be a character string with a formula as a function of
uppercase X, Y, and Z. The most common is the linear trend,
\code{"~ X + Y + Z"}. For ordinary kriging, use \code{"~1"}. If neither
\code{trend} nor \code{mean} are given, it is assumed that the global mean
is the mean of the data values.

If any point index is given in \code{force.interp}, the predicted mean
function will pass exactly through those points, but the predictive variance
will still be computed as usual. This is in contrast to what is usually
done by geostatistics softwares, which assign a variance of 0 to those
points.

\code{weights} can be used to "tune down" the effect of some data points.
The smaller the weight, the less effect a point will have on the predicted
funtion, to the limit that a weight of 0 filters the point completely.

Note that this implementation uses all the data provided to make predictions,
which may be too memory intensive for large datasets.
}
\seealso{
\code{\link{GP-class}}, \code{\link{SPGP-class}}
}

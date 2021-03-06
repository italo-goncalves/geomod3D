% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/lines3DDataFrame.R
\docType{methods}
\name{MergeSegments}
\alias{MergeSegments}
\alias{MergeSegments,lines3DDataFrame-method}
\title{MergeSegments}
\usage{
MergeSegments(object, ...)

\S4method{MergeSegments}{lines3DDataFrame}(object, by,
  keep = colnames(GetData(object)))
}
\arguments{
\item{by}{Column names to be used on merging.}

\item{keep}{Which columns to keep in the merged object.}
}
\value{
A \code{lines3DDataFrame} object.
}
\description{
Merging of redundant line segments to reduce the number of rows in a
\code{lines3DDataFrame} object.
}
\details{
This method merges contiguous line segments that are parallel and
for which all variables described in the \code{by} parameter are equal.
Variables described in the \code{keep} parameter will be averaged if numeric.
If any variable described in \code{keep} is a character of factor, new
columns will be created in the merged object that contain the relative
proportion of each unique element or factor in each merged segment.
}

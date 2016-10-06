## This demo script shows some of functionalities of the project and 
## a modeling example.
##
## Before running this script, make sure to download all source files and
## unzip the data file. The code assumes everything is in the same folder.
##
## It is recommended to run the script line by line in order to see what
## happens at each step.

#### Packages ####
library(rgl)
library(misc3d)
library(GA)
library(Rcpp)
library(Matrix)

#### Project's source files ####
sourceCpp("covC.cpp")
source('generics.R')
source('point3D.R')
source('line3D.R')
source('spatial3DDataFrame.R')
source('points3DDataFrame.R')
source('lines3DDataFrame.R')
source('grid3DDataFrame.R')
source('functions.R')
source('covarianceStructure3D.R')
source('geostat3D.R')
source('GP.R')
source('GP_LR.R')

#### Opening data ####
# set working directory to the source files location

# opening data
holes <- read.csv("data.csv", sep = ";", stringsAsFactors = F)

# a separate data frame is created which contains only the orientation data
dips <- holes[!is.na(holes$Dip),]

# calculation of strike (to be used later)
dips$Strike <- dips$Azimuth - 90
dips$Strike[dips$Strike < 0] <- dips$Strike[dips$Strike < 0] + 360

#### geomod3D objects ####
# Most of the project's objects consist of data frames with coordinates
# attached, similar to the objects from the sp package

# object with data on the geological classes
# "Horizon.up" and "Horizon.down" represent the layers on each side of
# a point. If the two values are equal, the point is inside a layer.
# If the values are different, the point is located in a boundary.
pointdata <- points3DDataFrame(holes[,2:4], holes[,-c(2,3,4)])
pointdata

# Object with orientation data
dipdata <- points3DDataFrame(dips[,2:4], dips[,-c(2,3,4)])
dipdata

# Object with unit vectors pointing along dip and strike. It uses dip and 
# strike angles to do the calculations. Note that each plane contains two
# vectors, so the output has twice the number of rows as the input.
dipvec <- getPlaneDirections(dipdata)
dipvec

#### Model training with a genetic algorithm ####
# A genetic algorithm is used to find the parameters with best the fit to 
# the data. The parameters are the covariance function's amplitude, range
# in the horizontal direction, range in the vertical direction (encoded as
# a multiple of the horizontal range), and the nugget or noise parameter.
#
# The covariance function is encoded in the "covarianceStructure3D" object.
# The model itself is the "GP_LR" object (Gaussian Process for Logistic
# Regression). It contains a separate Gaussian Process for each class and
# takes the covariance function and nugget as inputs, as well as the point
# and directional data. The fitness function to be maximized is the model's
# log-likelihood.
#
# This part can take some time.
fit <- function(x) {
        m <- covarianceStructure3D("gaussian", x[1], x[2], x[2], x[2]*x[3])
        gp <- GP_LR(pointdata, value1 = "Horizon.up", value2 = "Horizon.down",
                      model = m, strength = 1, nugget = x[4],
                      tangents = dipvec, reg.t = 1e-8)
        return(logLik(gp))
}
GA <- ga(type = "real-valued", fitness = fit,
         min = c(0.1, 100, 0.1, 0.05),
         max = c(5, 1000, 1, 5), pmutation = 0.3,
         popSize = 20, run = 10, monitor = T,
         names = c("amplitude", "range", "rmult","nugget"))
summary(GA)

#### Final model and prediction ####
# The model is built using the best fit parameters from the code above.
# For reproducibility, the parameters used in the paper are provided.
m <- covarianceStructure3D("gaussian", 0.65, 332, 332, 146)
m
gp <- GP_LR(pointdata, value1 = "Horizon.up", value2 = "Horizon.down", 
            strength = 1, model = m, nugget = 0.59, tangents = dipvec, 
            reg.t = 1e-8)
gp
logLik(gp)

# A regular, empty 3D grid is initialized
gr <- grid3DDataFrame(gridx = seq(-10100, -9700, 10), 
                      gridy = seq(-3900, -3100, 10),
                      gridz = seq(-1000, 100, 10), 
                      fields = "Horizon")
# It is then overwritten with the model's prediction (this can be slow)
gr <- predict(gp, gr, to = "Horizon", output.unc = T)
gr

#### Visualization of data and results ####
# Color pallette
hcol <- c(
        "purple3", "seagreen", "saddlebrown", "royalblue3", "orange2",
        "lightsalmon3", "lightgoldenrod3", "lawngreen", "lightblue", "darkorchid",
        "darkorange3", "forestgreen", "firebrick2", "dodgerblue3", "khaki2",
        "hotpink1", "honeydew2", "deepskyblue3", "gold", "cyan",
        "chartreuse", "burlywood3", "coral3", "aquamarine4", "maroon1",
        "olivedrab", "orangered", "slateblue1", "yellow3", "springgreen3",
        "turquoise4", "thistle"
)

# Visualization of data
drawPoints(pointdata, size = 10, by = "Horizon.up", 
           values = sort(unique(holes$Horizon.up)), col = hcol)
drawTangentPlanes(dipdata, size = 50)
axes3d()

# Visualization of model
# The countouring algorithm may produce some artifacts due to some layers'
# thickness being smaller than the grid spacing
for(i in 0:30){
        hor <- make3DArray(gr, paste0("Horizon..Horizon_", 
                                       sprintf("%02d", i),".ind"))
        contour3d(hor$value, level = 0, x = hor$x, y = hor$y, z = hor$z, 
                  color = hcol[i+1], alpha = 1, add=T)
}

# Visualization of uncertainty (relative predictive standard deviation)
# The "[[" extracts a single variable from the grid
hist(gr[["Horizon..uncertainty"]])

# Uncertainty contouring
unc <- make3DArray(gr, "Horizon..uncertainty")
contour3d(unc$value, level = c(0.5,0.75,0.99), 
          x = unc$x, y = unc$y, z = unc$z, 
          color = c("green", "yellow", "red"), alpha = 0.5, add=T)

# Model constrained by an uncertainty threshold
thr <- 0.75
mask <- unc$value <= thr
for(i in 0:30){
        hor <- make3DArray(gr, paste0("Horizon..Horizon_", 
                                      sprintf("%02d", i),".ind"))
        contour3d(hor$value, level = 0, x = hor$x, y = hor$y, z = hor$z, 
                  color = hcol[i+1], alpha = 1, add=T, mask = mask)
}
contour3d(unc$value, level = thr, 
          x = unc$x, y = unc$y, z = unc$z, 
          color = "gray50", alpha = 0.5, add=T)

# Reference surfaces
load("ref.Rdata")
for(i in 1:31){
        persp3d(sup[[i]]$x, sup[[i]]$y, sup[[i]]$value, add = T, col = hcol[i],
                xlim = c(-10100,-9700), ylim = c(-3800,-3200), 
                zlim = c(-1000, 100))
}
## This demo script shows some of functionalities of the project and 
## a modeling example.
##
## Before running this script, make sure to download all source files and
## unzip the data file. The code assumes everything is in the same folder.
##
## The sample data was exported from the SCAT module of software Move, which is 
## available for universities and research institutions through the Academic 
## Software Initiative (ASI). The data should be used only within this script.
##
## It is recommended to run the script line by line in order to see what
## happens at each step.

#### Packages ####
library(rgl)
library(misc3d)
library(Rcpp)
library(Matrix)
library(ggplot2)

#### Project's source files ####
sourceCpp("covC.cpp")
source('generics.R')
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
holes <- read.csv(unz("data.zip", "data.csv"), 
                  sep = ";", stringsAsFactors = F)

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

#### Model training ####
# The model's parameters (range and nugget effect) are found by maximizing
# the log-likelihood, as there are only two parameters, it suffices to do
# a grid search and plot the results. If one chooses to work with more 
# parameters (anisotropy, different parameters for each class, etc.), a
# more complex training method must be used.
#
# The covariance function is encoded in the "covarianceStructure3D" object.
# The model itself is the "GP_LR" object (Gaussian Process for Logistic
# Regression). It contains a separate Gaussian Process for each class and
# takes the covariance function and nugget as inputs, as well as the point
# and directional data. Valid covariances are the Gaussian and cubic functions.
# The Matérn class is also supported but did not gave good results with this 
# dataset.

# The following function plots the grid search results
plot_lik <- function(r, nug, lik){
        df <- data.frame(
                r = rep(r, times = length(nug)),
                nug = rep(nug, each = length(r)),
                lik = as.numeric(lik)
        )
        ggplot(df) + geom_tile(aes(x=r, y=nug, fill=lik)) + 
                scale_fill_gradientn(colours = rev(colorRamps::matlab.like2(20)),
                                     name = "Log-likelihood") +
                scale_x_continuous(labels = r, breaks = r, expand = c(0,0),
                                   name = "Range") + 
                scale_y_continuous(labels = nug, breaks = nug, expand = c(0,0)) + 
                ylab(expression(~sigma[0]^2))
}

# The amplitude parameter can be fixed
s12 <- 33^2/(16*32^2)

# grid search parameters
r <- seq(100,1000,100) # ranges to test
nug <- seq(0.02, 0.3, 0.02) # nugget values to test
lik <- matrix(nrow = length(r), ncol = length(nug)) # matrix to store the results

# covariance function (one of "gaussian", "cubic", "matern1", or "matern2")
covfun <- "gaussian"

# This part can take some time
for(i in seq_along(r)){
        for(j in seq_along(nug)){
                lik[i,j] <- logLik(GP_LR(pointdata, "Horizon.up", "Horizon.down", 
                                         model = covarianceStructure3D(
                                                 type = covfun, 
                                                 contribution = s12, 
                                                 maxrange = r[i]), 
                                         nugget = nug[j], 
                                         tangents = dipvec, 
                                         reg.t = 1e-8))
        }
}
plot_lik(r, nug, lik)


#### Final model and prediction ####
# The model is built using a combination of parameters from the plot above.
# For reproducibility, the parameters used in the paper are provided.
# (reg.t is the regularization parameter for the tangent data)
m <- covarianceStructure3D("gaussian", s12, 900)
m
gp <- GP_LR(pointdata, value1 = "Horizon.up", value2 = "Horizon.down", 
            model = m, nugget = 0.22, tangents = dipvec, reg.t = 1e-8)
gp


# A regular, empty 3D grid is initialized
gr <- grid3DDataFrame(gridx = seq(-10100, -9700, 10), 
                      gridy = seq(-3900, -3100, 10),
                      gridz = seq(-1000, 100, 10), 
                      fields = "Horizon")
# It is then overwritten with the model's prediction (this can be slow)
gr <- predict(gp, gr, to = "Horizon")
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


# Reference surfaces
load(unz("data.zip", "ref.Rdata"))
for(i in 1:31){
        persp3d(sup[[i]]$x, sup[[i]]$y, sup[[i]]$value, add = T, col = hcol[i],
                xlim = c(-10100,-9700), ylim = c(-3800,-3200), 
                zlim = c(-1000, 100))
}
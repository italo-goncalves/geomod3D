#### points3DDataFrame class ####
points3DDataFrame <- setClass(
        "points3DDataFrame",
        contains = "spatial3DDataFrame",
        validity = function(object) {
                if (length(object@coords) !=
                    nrow(object@data)) {
                        stop(
                                paste(
                                        "Number of point3D objects does not",
                                        "match number of observations"
                                )
                        )
                }
                if (!all(rapply(object@coords,class) == "point3D")) {
                        stop(
                                paste(
                                        "Invalid object in",
                                        "points list. All",
                                        "objects must be of",
                                        "class 'point3D'"
                                )
                        )
                }
        }
)

#### initialization ####
setMethod(
        f = "initialize",
        signature = "points3DDataFrame",
        definition = function(.Object, coords, df){
                if(class(coords) %in% c("matrix", "Matrix", "data.frame", 
                                        "tbl_df")) {
                        if(!(ncol(coords) %in% c(2,3)))
                                stop("Invalid number of dimensions")
                        coords <- apply(coords, 1, function(x) point3D(x))
                        .Object@coords <- coords
                }else if(class(coords) == "list"){
                        .Object@coords <- coords
                }else{
                        stop("Invalid format for coordinates")
                }
                .Object@data <- df
                # bounding box
                if(nrow(df) > 0){
                        points_df <- getCoords(.Object, "data.frame")
                        bbox <- as.matrix(rbind(
                                apply(points_df,2,min),
                                apply(points_df,2,max)))
                        rownames(bbox) <- c("min","max")
                        .Object@bbox <- bbox
                }else .Object@bbox <- matrix(0, 2, 3)
                # end
                validObject(.Object)
                return(.Object)
        }
)

#### getCoords ####
setMethod(
        f = "getCoords",
        signature = "points3DDataFrame",
        definition = function(object, as = "list"){
                if(as == "list"){
                        return(object@coords)
                }else if(as == "data.frame"){
                        points_list <- object@coords
                        df <- matrix(nrow = length(points_list), ncol = 3)
                        for(i in seq(length(points_list))){
                                df[i,] <- getCoords(points_list[[i]])
                        }
                        df <- as.data.frame(df)
                        colnames(df) <- c("X","Y","Z")
                        return(df)
                }else if(as == "matrix"){
                        points_list <- object@coords
                        df <- matrix(nrow = length(points_list), ncol = 3)
                        for(i in seq(length(points_list))){
                                df[i,] <- getCoords(points_list[[i]])
                        }
                        colnames(df) <- c("X","Y","Z")
                        return(df)
                }else{
                        stop("Invalid value for 'as' argument")
                }
        }
)

#### setAs ####
setAs("NULL", "points3DDataFrame", function(from, to)
        new(to, list(), data.frame()))


#### rbind, cbind equivalents ####
setMethod("bindPoints", c("points3DDataFrame","points3DDataFrame"),
          function(x, y){
                  points_list <- c(getCoords(x), getCoords(y))
                  df <- merge(getData(x), getData(y), all=T, sort=F)
                  return(points3DDataFrame(points_list,df))
          })

#### reScale ####
setMethod(
        f = "reScale",
        signature = "points3DDataFrame",
        definition = function(object, 
                              old_range = boundingBox(object), 
                              new_range = matrix(rep(c(0,1),3),2,3)){
                points_list <- getCoords(object)
                df <- getData(object)
                points_list <- lapply(
                        points_list,
                        function(x) reScale(x, old_range, new_range)
                )
                return(points3DDataFrame(points_list,df))
        }
)

setMethod(
        f = "pointify",
        signature = "points3DDataFrame",
        definition = function(object) object
)

#### drawPoints ####
setMethod(
        f = "drawPoints",
        signature = "points3DDataFrame",
        definition = function(object, by, values, col, size, alpha = 1){
                # pacakges
                require(rgl)
                # require(grDevices)
                # require(plotrix)
                
                # setup
                coords <- getCoords(object, "matrix")
                df <- getData(object)
                N <- nrow(object)
                objval <- df[,by]
                if(length(size) < N) size <- rep(size, length.out = N)
                if(length(alpha) < N) alpha <- rep(alpha, length.out = N)
                
                # pallette
                if(class(objval) == "numeric"){ # continuous variable
                        colorsc <- find_color_cont(objval, 
                                                   rng = range(values),
                                                   col = col)
                }else{ # categorical variable
                        names(col) <- values
                        colorsc <- col[objval]
                }
                
                # plotting
                spheres3d(coords, radius = size/2,
                          color = colorsc, alpha = alpha)
        }
)

#### drawTangentPlanes ####
setMethod(
        f = "drawTangentPlanes",
        signature = "points3DDataFrame",
        definition = function(object, size, dip = "Dip", strike = "Strike",
                              col = "yellow"){
                require(rgl)
                # setup
                N <- nrow(object)
                normalvec <- getData(getNormals(object, dip, strike))
                coords <- getCoords(object, "matrix")
                if(length(col) < N) col <- rep(col, length.out = N)
                
                # drawing planes
                for(i in seq(nrow(object))){
                        cylcoords <- rbind(coords[i,] + normalvec[i,] * size/1000,
                                           coords[i,] - normalvec[i,] * size/1000)
                        shade3d(
                                cylinder3d(cylcoords, radius = size/2, 
                                           sides = 128, closed = -2),
                                col = col[i]
                        )
                }
        }
)

#### drawTangentLines ####
setMethod(
        f = "drawTangentLines",
        signature = "points3DDataFrame",
        definition = function(object, size, dX = "dX", dY = "dY", dZ = "dZ",
                              col = "yellow"){
                require(rgl)
                # setup
                N <- nrow(object)
                coords <- getCoords(object, "matrix")
                dirs <- getData(object[c(dX, dY, dZ)])
                if(length(col) < N) col <- rep(col, length.out = N)
                
                # drawing lines
                for(i in seq(N)){
                        cylcoords <- rbind(coords[i,] + dirs[i,] * size/2,
                                           coords[i,] - dirs[i,] * size/2)
                        shade3d(
                                cylinder3d(cylcoords, radius = size/10, 
                                           sides = 32, closed = -2),
                                col = col[i]
                        )
                }
        }
)

#### contactsBestFitPlane ####
setMethod(
        f = "contactsBestFitPlane",
        signature = "points3DDataFrame",
        definition = function(object, up, down, clusters, normalize = F){
                # Fernández, 2005
                # setup
                df <- getData(object)
                
                # finding contacts
                is.contact <- df[, up] != df[, down]
                coords <- getCoords(object, "matrix")[is.contact, ]
                
                # clustering
                d <- dist(coords)
                hc <- hclust(d, "single")
                cl <- cutree(hc, k = clusters)
                
                # moment of inertia
                vec1 <- vec2 <- as.data.frame(matrix(0, nrow(coords), 6))
                colnames(vec1) <- colnames(vec2) <- 
                        c("dX", "dY", "dZ", "M", "K", "cluster")
                keep <- logical(nrow(coords))
                for(i in seq(clusters)){
                        id <- which(cl == i)
                        if(length(id) >= 3){
                                tmpcoords <- coords[id, , drop = F]
                                tmpcoords <- apply(tmpcoords, 2, 
                                                   function(cl) cl - mean(cl))
                                if(normalize) 
                                        tmpcoords <- 
                                        t(apply(tmpcoords, 1, 
                                                function(rw) rw/sqrt(sum(rw^2))))
                                s <- svd(tmpcoords)
                                vec1[id,1:3] <- matrix(s$v[,1], 
                                                       length(id), 3, byrow = T)
                                vec2[id,1:3] <- matrix(s$v[,2], 
                                                       length(id), 3, byrow = T)
                                d <- s$d^2
                                vec1[id, "K"] <- vec2[id, "K"] <- log(d[1]/d[3])
                                vec1[id, "M"] <- vec2[id, "M"] <- 
                                        log(d[1]/d[2]) / log(d[2]/d[3])
                                vec1[id, "cluster"] <- vec2[id, "cluster"] <- i
                                keep[id] <- TRUE
                        }
                }
                
                # result
                return(points3DDataFrame(
                        rbind(coords[keep,], coords[keep,]),
                        rbind(vec1[keep,], vec2[keep,])
                ))
                
        }
)

#### classify ####
setMethod(
        f = "classify",
        signature = c("points3DDataFrame", "points3DDataFrame"),
        definition = function(x, y, value1, value2 = value1, to = value1, 
                              strength, model, nugget, verbose = T, 
                              output.unc = F,
                              tangents = NULL, dip = "Dip", strike = "Strike"){
                # setup
                Ndata <- nrow(x)
                xdata <- getData(x)
                ydata <- getData(y)
                categories <- unique(c(xdata[,value1], xdata[,value2]))
                ncat <- length(categories)
                # y[".mean"] <- -1
                if(!is.null(tangents)) 
                        tangents <- getPlaneDirections(tangents, dip, strike)

                # indicator matrix
                indmat <- matrix(NA, nrow(y), ncat)
                varmat <- matrix(NA, nrow(y), ncat)

                # indicator kriging
                for(item in categories){
                        # indicator
                        ind <- matrix(-strength, nrow(x), 2) # negative
                        ind[xdata[,value1] == item, 1] <- strength # positive 1
                        ind[xdata[,value2] == item, 2] <- strength # positive 2
                        x[".ind"] <- rowMeans(ind) # contacts get 0
                        x[".nugget"] <- nugget
                        x[getData(x)$.ind == 0, ".nugget"] <- 1e-9 # regularization
                        # kriging
                        if(verbose){
                                cat("Classifying ", item, " ... \n",
                                    sep = "")
                        }
                        gp <- GP(x, model, value = ".ind", nugget = ".nugget",
                                 mean = -1, tangents = tangents)
                        if(output.unc){
                                tempgr <- predict(gp, y, output.var = T)
                                indmat[, which(categories == item)] <-
                                        unlist(getData(tempgr[".ind"]))
                                varmat[, which(categories == item)] <-
                                        unlist(getData(tempgr[".ind.var"]))
                        }
                        else{
                                tempgr <- predict(gp, y, output.var = F)
                                indmat[, which(categories == item)] <-
                                        unlist(getData(tempgr[".ind"]))
                        }
                }


                # classification (one vs all)
                # code_matrix <- matrix(-1,ncat,ncat) + diag(2,ncat,ncat)
                # indmat <- indmat %*% code_matrix
                calcprob <- function(rw){
                        # preventing NaNs
                        rw[rw > 20] <- 20
                        rw[rw < -20] <- -20
                        # randomly breaking ties
                        rw <- rw + rnorm(length(rw), sd = 1e-6)
                        # probabilities
                        exp(rw) / sum(exp(rw))
                }
                probmat <- t(apply(indmat, 1, calcprob))
                ids <- apply(probmat, 1, which.max)
                ydata[,to] <- categories[ids]

                # probabilities
                colnames(probmat) <- paste0(to,"..",
                                            make.names(categories),
                                            ".prob")
                ydata[,colnames(probmat)] <- as.data.frame(probmat)

                # indicators
                # distance from border calculated from probabilities
                indmat <- t(apply(probmat, 1, function(x){
                        y <- log(x)
                        ysort <- sort(y, decreasing = T)
                        y - mean(ysort[1:2])
                }))
                inddf <- data.frame(indmat)
                colnames(inddf) <- paste0(to,"..",
                                          make.names(categories),".ind")
                ydata[,colnames(inddf)] <- inddf

                # uncertainty (geometric mean of entropy and relative variance)
                if(output.unc){
                        if(class(model) != "list") model <- list(model)
                        totvar <- sum(sapply(model, function(m) m@contribution))
                        entropy <- rowSums(- probmat * log(probmat)) / log(ncat)
                        rel_var <- rowMeans(varmat) / totvar
                        u <- sqrt(entropy * rel_var)
                        ydata[, paste0(to,"..uncertainty")] <- u
                }

                # end
                y@data <- ydata
                return(y)
        }
)

setMethod(
        f = "classify",
        signature = c("points3DDataFrame", "missing"),
        definition = function(x, y, value1, value2 = value1, 
                              model, strength, nugget,
                              verbose = T,
                              tangents = NULL, dip = "Dip", strike = "Strike",
                              parts = 5, seed = NULL){
                # setup
                xdata <- getData(x)
                categories <- unique(c(xdata[,value1], xdata[,value2]))
                ncat <- length(categories)
                id_points <- which(xdata[,value1] == xdata[,value2])
                if(!is.null(seed)) set.seed(seed)
                CVpart <- sample.int(parts, size = nrow(x),
                                     replace = T)
                pred <- character(length = nrow(x))

                # cross validation
                for(i in seq(parts)){
                        if(verbose) cat("Cross validation partition", i,
                                        "of", parts, "\n")
                        # xids <- id_points[CVpart != i]
                        # yids <- id_points[CVpart == i]
                        # xtemp <- bindPoints(x[-id_points,], x[xids,])
                        # ytemp <- x[yids, ]
                        xids <- CVpart != i
                        yids <- CVpart == i
                        ytemp <- classify(x[xids,], x[yids,], value1 = value1,
                                          value2 = value2, model = model,
                                          strength = strength, nugget = nugget,
                                          verbose = verbose,
                                          to = ".pred",
                                          tangents = tangents, dip = dip,
                                          strike = strike)
                        pred[yids] <- unlist(getData(ytemp[".pred"]))
                }
                cv <- mean(pred[id_points] == xdata[id_points, value1])
                return(cv)
        }
)

# setMethod(
#         f = "classify_ps",
#         signature = c("points3DDataFrame", "points3DDataFrame"),
#         definition = function(x, y, value1, value2 = value1, to = value1, model,
#                               nugget, verbose = T, output.unc = F,
#                               tangents = NULL, dip = "Dip", strike = "Strike",
#                               pseudo){
#                 # setup
#                 Ndata <- nrow(x)
#                 xdata <- getData(x)
#                 ydata <- getData(y)
#                 categories <- unique(c(xdata[,value1], xdata[,value2]))
#                 ncat <- length(categories)
#                 # y[".mean"] <- -1
#                 if(!is.null(tangents)) 
#                         tangents <- getPlaneDirections(tangents, dip, strike)
#                 
#                 # indicator matrix
#                 indmat <- matrix(NA, nrow(y), ncat)
#                 varmat <- matrix(NA, nrow(y), ncat)
#                 
#                 # indicator kriging
#                 for(item in categories){
#                         # indicator
#                         ind <- matrix(-1, nrow(x), 2) # negative
#                         ind[xdata[,value1] == item, 1] <- 1 # positive 1
#                         ind[xdata[,value2] == item, 2] <- 1 # positive 2
#                         x[".ind"] <- rowMeans(ind) # contacts get 0
#                         x[".nugget"] <- nugget
#                         x[getData(x)$.ind == 0, ".nugget"] <- 1e-9 # regularization
#                         # kriging
#                         if(verbose){
#                                 cat("Classifying ", item, " ... \n",
#                                     sep = "")
#                         }
#                         gp <- SPGP(x, model, value = ".ind", nugget = ".nugget",
#                                  mean = -1, 
#                                  pseudo_inputs = pseudo[,1:3],
#                                  pseudo_nugget = pseudo[,4])
#                         if(output.unc){
#                                 tempgr <- predict(gp, y, output.var = T)
#                                 indmat[, which(categories == item)] <-
#                                         unlist(getData(tempgr[".ind"]))
#                                 varmat[, which(categories == item)] <-
#                                         unlist(getData(tempgr[".ind.var"]))
#                         }
#                         else{
#                                 tempgr <- predict(gp, y, output.var = F)
#                                 indmat[, which(categories == item)] <-
#                                         unlist(getData(tempgr[".ind"]))
#                         }
#                 }
#                 
#                 
#                 # classification (one vs all)
#                 code_matrix <- matrix(-1,ncat,ncat) + diag(2,ncat,ncat)
#                 indmat <- indmat %*% code_matrix
#                 calcprob <- function(rw){
#                         # preventing NaNs
#                         rw[rw > 10] <- 10
#                         rw[rw < -10] <- -10
#                         # randomly breaking ties
#                         rw <- rw + rnorm(length(rw), sd = 1e-6)
#                         # probabilities
#                         exp(rw) / sum(exp(rw))
#                 }
#                 probmat <- t(apply(indmat, 1, calcprob))
#                 ids <- apply(probmat, 1, which.max)
#                 ydata[,to] <- categories[ids]
#                 
#                 # probabilities
#                 colnames(probmat) <- paste0(to,"..",
#                                             make.names(categories),
#                                             ".prob")
#                 ydata[,colnames(probmat)] <- as.data.frame(probmat)
#                 
#                 # indicators
#                 # distance from border calculated from probabilities
#                 indmat <- t(apply(probmat, 1, function(x){
#                         y <- log(x)
#                         ysort <- sort(y, decreasing = T)
#                         y - mean(ysort[1:2])
#                 }))
#                 inddf <- data.frame(indmat)
#                 colnames(inddf) <- paste0(to,"..",
#                                           make.names(categories),".ind")
#                 ydata[,colnames(inddf)] <- inddf
#                 
#                 # uncertainty (geometric mean of entropy and relative variance)
#                 if(output.unc){
#                         if(class(model) != "list") model <- list(model)
#                         totvar <- sum(sapply(model, function(m) m@contribution))
#                         entropy <- rowSums(- probmat * log(probmat)) / log(ncat)
#                         rel_var <- rowMeans(varmat) / totvar
#                         u <- sqrt(entropy * rel_var)
#                         ydata[, paste0(to,"..uncertainty")] <- u
#                 }
#                 
#                 # end
#                 y@data <- ydata
#                 return(y)
#         }
# )
# 
# setMethod(
#         f = "classify_ps",
#         signature = c("points3DDataFrame", "missing"),
#         definition = function(x, y, value1, value2 = value1, model, nugget,
#                               verbose = T, pseudo,
#                               tangents = NULL, dip = "Dip", strike = "Strike",
#                               parts = 5, seed = NULL){
#                 # setup
#                 xdata <- getData(x)
#                 categories <- unique(c(xdata[,value1], xdata[,value2]))
#                 ncat <- length(categories)
#                 id_points <- which(xdata[,value1] == xdata[,value2])
#                 if(!is.null(seed)) set.seed(seed)
#                 CVpart <- sample.int(parts, size = length(id_points),
#                                      replace = T)
#                 pred <- character(length = length(id_points))
#                 
#                 # cross validation
#                 for(i in seq(parts)){
#                         if(verbose) cat("Cross validation partition", i,
#                                         "of", parts, "\n")
#                         xids <- id_points[CVpart != i]
#                         yids <- id_points[CVpart == i]
#                         xtemp <- bindPoints(x[-id_points,], x[xids,])
#                         ytemp <- x[yids, ]
#                         ytemp <- classify_ps(xtemp, ytemp, value1 = value1,
#                                           value2 = value2, model = model,
#                                           nugget = nugget,
#                                           verbose = verbose,
#                                           to = ".pred",
#                                           pseudo = pseudo)
#                         pred[CVpart == i] <- unlist(getData(ytemp[".pred"]))
#                 }
#                 cv <- mean(pred == xdata[id_points, value1])
#                 return(cv)
#         }
# )

# # compositional
# setMethod(
#         f = "classify",
#         signature = c("points3DDataFrame", "points3DDataFrame"),
#         definition = function(x, y, value1, value2 = value1, to = value1, model, 
#                               nugget, verbose = T, prob.smooth = 0.01,
#                               tangents = NULL, dip = "Dip", strike = "Strike",
#                               output.unc = F){
#                 # setup
#                 Ndata <- nrow(x)
#                 xdata <- getData(x)
#                 ydata <- getData(y)
#                 categories <- unique(c(xdata[,value1], xdata[,value2]))
#                 ncat <- length(categories)
#                 y[".mean"] <- -1
#                 if(!is.null(tangents)) 
#                         tangents <- getPlaneDirections(tangents, dip, strike)
#                 
#                 # indicator matrix
#                 indmat <- matrix(0, Ndata, ncat)
#                 for(i in seq_along(categories)){
#                         id <- which(
#                                 xdata[,value1] == categories[i] |
#                                         xdata[,value2] == categories[i]
#                         )
#                         indmat[id,i] <- 1 - prob.smooth
#                 }
#                 indmat <- t(apply(indmat, 1, function(rw){
#                         id <- which(rw > 0)
#                         rw[id] <- rw[id] / length(id)
#                         rw[-id] <- prob.smooth / (ncat - length(id))
#                         log(rw) - mean(log(rw))
#                 }))
#                 indpred <- matrix(0, nrow(y), ncat)
#                 
#                 # nugget
#                 nugget <- rep(nugget, Ndata)
#                 id <- which(xdata[,value1] != xdata[,value1])
#                 nugget[id] <- 1e-9
#                 
#                 # pre-processing
#                 x[".ind"] <- indmat[,1]
#                 gp <- GP(x, model, value = ".ind", nugget = nugget,
#                          tangents = tangents)
#                 Ntang <- ifelse(is.null(tangents), 0, nrow(tangents))
#                 
#                 # Gaussian Process regression
#                 for(i in seq_along(categories)){
#                         if(verbose){
#                                 cat("Classifying ", categories[i], " ... \n", 
#                                     sep = "")
#                         }
#                         if(i > 1){
#                                 # changing variable
#                                 tmp <- c(indmat[,i], rep(0, Ntang))
#                                 gp@w_value <- apply(gp@w_var, 1, function(rw){
#                                         sum(rw * tmp)
#                                 })
#                         }
#                         if(i == ncat & output.unc){
#                                 tempgr <- predict(gp, y, output.var = T)
#                                 tempdata <- getData(tempgr)
#                                 indpred[,i] <- tempdata[, ".ind"]
#                                 var_pred <- tempdata[, ".ind.var"]
#                         }
#                         else{
#                                 tempgr <- predict(gp, y, output.var = F)
#                                 tempdata <- getData(tempgr)
#                                 indpred[,i] <- tempdata[, ".ind"]
#                         }
#                 }
#                 
#                 
#                 # classification
#                 probmat <- t(apply(indpred, 1, function(rw) exp(rw)/sum(exp(rw))))
#                 ids <- apply(probmat, 1, which.max)
#                 ydata[,to] <- categories[ids]
#                 
#                 # probabilities
#                 colnames(probmat) <- paste0(to,"..",
#                                             make.names(categories),
#                                             ".prob")
#                 ydata[,colnames(probmat)] <- as.data.frame(probmat)
#                 
#                 # indicators
#                 # distance from border calculated from probabilities
#                 indpred <- t(apply(probmat, 1, function(x){
#                         y <- log(x)
#                         ysort <- sort(y, decreasing = T)
#                         y - mean(ysort[1:2])
#                 }))
#                 inddf <- data.frame(indpred)
#                 colnames(inddf) <- paste0(to,"..",
#                                           make.names(categories),".ind")
#                 ydata[,colnames(inddf)] <- inddf
#                 
#                 # uncertainty (geometric mean of entropy and relative variance)
#                 if(output.unc){
#                         if(class(model) != "list") model <- list(model)
#                         totvar <- sum(sapply(model, function(m) m@contribution))
#                         entropy <- rowSums(- probmat * log(probmat)) / log(ncat)
#                         var_pred <- var_pred / totvar
#                         u <- sqrt(entropy * var_pred)
#                         ydata[, paste0(to,"..uncertainty")] <- u
#                 }
#                 
#                 # end
#                 y@data <- ydata
#                 return(y)
#         }
# )
# 
# setMethod(
#         f = "classify",
#         signature = c("points3DDataFrame", "missing"),
#         definition = function(x, y, value1, value2 = value1, model, nugget,
#                               verbose = T, prob.smooth = 0.01,
#                               tangents = NULL, dip = "Dip", strike = "Strike",
#                               parts = 5, seed = NULL){
#                 # setup
#                 xdata <- getData(x)
#                 categories <- unique(c(xdata[,value1], xdata[,value2]))
#                 ncat <- length(categories)
#                 id_points <- which(xdata[,value1] == xdata[,value2])
#                 if(!is.null(seed)) set.seed(seed)
#                 CVpart <- sample.int(parts, size = length(id_points),
#                                      replace = T)
#                 pred <- character(length = length(id_points))
#                 
#                 # cross validation
#                 for(i in seq(parts)){
#                         if(verbose) cat("Cross validation partition", i,
#                                         "of", parts, "\n")
#                         xids <- id_points[CVpart != i]
#                         yids <- id_points[CVpart == i]
#                         xtemp <- bindPoints(x[-id_points,], x[xids,])
#                         ytemp <- x[yids, ]
#                         ytemp <- classify(xtemp, ytemp, value1 = value1,
#                                           value2 = value2, model = model,
#                                           nugget = nugget,
#                                           verbose = verbose,
#                                           to = ".pred", 
#                                           prob.smooth = prob.smooth,
#                                           tangents = tangents, dip = dip, 
#                                           strike = strike)
#                         pred[CVpart == i] <- unlist(getData(ytemp[".pred"]))
#                 }
#                 cv <- mean(pred == xdata[id_points, value1])
#                 return(cv)
#         }
# )
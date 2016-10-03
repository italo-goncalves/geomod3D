

#### elipsóide ####
asymmetry3D <- function(maxrange, midrange, minrange = midrange,
                        azimuth = 0, dip = 0, rake = 0, radians = F){
        # conversion to radians
        if(!radians){
                azimuth <- azimuth * pi/180
                dip <- dip * pi/180
                rake <- rake * pi/180
        }
        
        # conversion to mathematical coordinates
        dip <- -dip
        r <- c(midrange, maxrange, minrange)
        
        # rotation matrix
        Rx <- diag(1,3,3)
        Rx[c(1,3),c(1,3)] <- c(cos(rake),-sin(rake),sin(rake),cos(rake))
        Ry <- diag(1,3,3)
        Ry[2:3,2:3] <- c(cos(dip),sin(dip),-sin(dip),cos(dip))
        Rz <- diag(1,3,3)
        Rz[1:2,1:2] <- c(cos(azimuth),-sin(azimuth),
                         sin(azimuth),cos(azimuth))
        A <- diag(r,3,3)
        
        return(Rz %*% Ry %*% Rx %*% A)
}

#### covariância ####
setMethod(
        f = "covariance_matrix",
        signature = c("lines3DDataFrame", "lines3DDataFrame"),
        definition = function(x, y, model, value1, value2 = value1,
                              covariance = F){
                # covariance model
                if(length(model) == 1 & class(model) != "list"){
                        model <- list(model)
                }
                if(!all(rapply(model,class) == "covarianceStructure3D")){
                        stop("model must be of class 'covarianceStructure3D'")
                }
                
                # cleaning NAs
                na1 <- is.na(getData(x[value1]))
                x1 <- x[!na1,value1]
                Nrows <- nrow(x1)
                na2 <- is.na(getData(x[value2]))
                x2 <- x[!na2,value2]
                Ncols <- nrow(x2)
                
                # line discretization
                cat("Discretizing lines...\n")
                parts <- 10
                x1[".id"] <- seq(Nrows)
                x2[".id"] <- seq(Ncols)
                x1p <- pointify(x1, seq(0,1,1/parts))
                x2p <- pointify(x2, seq(0,1,1/parts))
                points1 <- getCoords(x1p,"matrix")
                points2 <- getCoords(x2p,"matrix")
                ids1 <- unlist(getData(x1p[".id"]))
                ids2 <- unlist(getData(x2p[".id"]))
                
                # covariance matrix
                cat("Building covariance matrix...\n")
                Kfull <- matrix(0,Nrows*(parts+1),Ncols*(parts+1))
                for(md in model){
                        Kfull <- Kfull + variogram3D(points1,
                                                     points2, 
                                                     model = md, 
                                             covariance = covariance)
                }
                K <- matrix(0,Nrows,Ncols)
                for(i in seq(Nrows)){
                        for(j in seq(Ncols)){
                                K[i,j] <- mean(Kfull[ids1==i,ids2==j])
                        }
                }

                return(K)
        }
)

setMethod(
        f = "covariance_matrix",
        signature = c("lines3DDataFrame","points3DDataFrame"),
        definition = function(x, y, model, covariance = F){
                # covariance model
                if(length(model) == 1 & class(model) != "list"){
                        model <- list(model)
                }
                if(!all(rapply(model,class) == "covarianceStructure3D")){
                        stop("model must be of class 'covarianceStructure3D'")
                }
                
                # cleaning NAs
                # na1 <- is.na(getData(x[value1]))
                # x1 <- x[!na1,value1]
                Nrows <- nrow(x)
                Ncols <- nrow(y)
                
                # line discretization
                cat("Discretizing lines...\n")
                parts <- 10
                x1[".id"] <- seq(Nrows)
                x1p <- pointify(x1, seq(0,1,1/parts))
                points1 <- getCoords(x1p,"matrix")
                ids1 <- unlist(getData(x1p[".id"]))
                
                # covariance matrix
                cat("Building covariance matrix...\n")
                points2 <- getCoords(y, "matrix")
                Kfull <- matrix(0,Nrows*(parts+1),Ncols)
                for(md in model){
                        Kfull <- Kfull + variogram3D(points1,
                                                     points2, 
                                                     model = md, 
                                                     covariance = covariance)
                }
                K <- matrix(0,Nrows,Ncols)
                for(i in seq(Nrows)){
                        for(j in seq(Ncols)){
                                K[i,j] <- mean(Kfull[ids1==i,j])
                        }
                }
                
                return(K)
        }
)

setMethod(
        f = "covariance_matrix",
        signature = c("points3DDataFrame","points3DDataFrame"),
        definition = function(x, y, model, covariance = F){
                # setup
                Nrows <- nrow(x)
                Ncols <- nrow(y)
                
                # covariance model
                if(length(model) == 1 & class(model) != "list"){
                        model <- list(model)
                }
                if(!all(rapply(model,class) == "covarianceStructure3D")){
                        stop("model must be of class 'covarianceStructure3D'")
                }

                # covariance matrix
                u <- getCoords(x, as = "matrix")
                v <- getCoords(y, as = "matrix")
                K <- matrix(0,Nrows,Ncols)
                if(covariance) for(md in model) K <- K + md@covfun(u, v)
                else for(md in model) K <- K + md@varfun(u, v)
                return(K)
                
        }
)

setMethod(
        f = "covariance_matrix_d1", # value/derivative covariance
        signature = c("points3DDataFrame", "points3DDataFrame"),
        definition = function(x, tangents, model, covariance = F){
                # setup
                Ndata <- nrow(x)
                Ntang <- nrow(tangents)
                xcoords <- getCoords(x, "matrix")
                tcoords <- getCoords(tangents, "matrix")

                # covariance model
                if(length(model) == 1 & class(model) != "list"){
                        model <- list(model)
                }
                if(!all(rapply(model,class) == "covarianceStructure3D")){
                        stop("model must be of class 'covarianceStructure3D'")
                }

                # dip and strike vectors
                vec <- as.matrix(getData(tangents[c("dX","dY","dZ")]))

                # covariance matrix
                K <- matrix(0, Ndata, Ntang)
                if(covariance) for(md in model) 
                        K <- K + md@covd1(xcoords, tcoords, vec)
                else for(md in model) 
                        K <- K + md@vard1(xcoords, tcoords, vec)
                return(K)
                
        }
)

setMethod(
        f = "covariance_matrix_d2", # derivative/derivative covariance
        signature = "points3DDataFrame",
        definition = function(tangents, model, covariance = F){
                # setup
                Ntang <- nrow(tangents)
                tcoords <- getCoords(tangents, "matrix")
                
                # covariance model
                if(length(model) == 1 & class(model) != "list"){
                        model <- list(model)
                }
                if(!all(rapply(model,class) == "covarianceStructure3D")){
                        stop("model must be of class 'covarianceStructure3D'")
                }
                
                # tangent vectors
                dirvecs <- as.matrix(getData(tangents[c("dX","dY","dZ")]))

                # covariance matrix
                K <- matrix(0, Ntang, Ntang)
                if(covariance) for(md in model) {
                        K <- K + md@covd2(tcoords, tcoords, dirvecs, dirvecs)
                }
                else for(md in model) {
                        K <- K + md@vard2(tcoords, tcoords, dirvecs, dirvecs)
                }
                return(K)
                
        }
)

## trend matrix
setMethod(
        f = "trend_matrix",
        signature = c("points3DDataFrame", "character"),
        definition = function(x, trend){
                trend <- as.formula(trend)
                TR <- model.matrix(trend, getCoords(x, "data.frame"))
                return(TR)
        }
)

setMethod(
        f = "trend_matrix_d1",
        signature = c("points3DDataFrame", "character"),
        definition = function(x, trend){
                # setup
                coords <- getCoords(x, "data.frame")
                trend <- as.formula(trend)
                # trend gradient
                dTX <- model.matrix(deriv_formula(trend, "X"), coords)
                dTY <- model.matrix(deriv_formula(trend, "Y"), coords)
                dTZ <- model.matrix(deriv_formula(trend, "Z"), coords)
                dTX <- rbind(dTX, dTX)
                dTY <- rbind(dTY, dTY)
                dTZ <- rbind(dTZ, dTZ)
                # dip and strike vectors
                vec <- as.matrix(getData(x[c("dX","dY","dZ")]))
                # directional derivatives
                dT <- dTX
                for(i in seq(dim(dTX)[2])){
                        dT[,i] <- rowSums(
                                cbind(dTX[,i], dTY[,i], dTZ[,i]) * vec
                        )
                }
                return(dT)
        }
)

#### kriging ####
setMethod(
        f = "krig3D",
        signature = c("ANY","ANY"),
        definition = function(x, y, model, value, to = value, 
                              nugget = 0, mean = NULL,
                              trend = NULL, weights = F,
                              tangents = NULL, dip = "Dip", strike = "Strike"){
                require(Matrix)
                
                ## setup
                # verifying value
                xdata <- getData(x)
                if(!(value %in% colnames(xdata))){
                        stop(paste0("Inexistent value '",value,"'\n"))
                }
                
                # filtering NAs 
                na_ids <- !is.na(getData(x)[,value]) 
                
                # counting dimensions
                Ndata <- nrow(x[na_ids,])
                Ngrid <- nrow(y)
                Ntang <- ifelse(is.null(tangents), 0, nrow(tangents))
                Ntrend <- 0
                
                # parsing nugget
                if(length(nugget) == 1){
                        if(class(nugget) == "character"){
                                nugget <- xdata[,nugget]
                                nugget <- diag(nugget[na_ids])
                        }else{
                                nugget <- diag(nugget,Ndata,Ndata)
                        }
                }else if(is.null(dim(nugget)) || min(dim(nugget)) == 1){
                        nugget <- diag(nugget[na_ids])
                }
                
                # parsing covariance model
                if(class(model) != "list") model <- list(model)
                if(!all(rapply(model,class) == "covarianceStructure3D")){
                        stop("model must be of class 'covarianceStructure3D'")
                }
                totvar <- sum(sapply(model, function(m) m@contribution))
                covariance <- is.null(trend)
                
                # slicing y to save memory
                maxgrid <- 1000 # optimize this
                Nslice <- ceiling(Ngrid / maxgrid)
                
                ## building system of equations
                # covariance matrix
                K <- covariance_matrix(x, x, model = model, 
                                       covariance = covariance)
                
                
                # trend matrix
                if(is.null(trend)){
                        TR <- matrix(0,Ndata,0)
                        K <- K + nugget
                }else{
                        TR <- trend_matrix(x, trend)
                        K <- K - nugget
                        Ntrend <- dim(TR)[2]
                }
                
                # tangents
                dK <- matrix(0, Ndata, 0)
                ddK <- matrix(0, 0, 0)
                # dK0 <- matrix(0, 0, Ngrid) ###
                dT <- matrix(0, 2*Ntang, Ntrend)
                if(!is.null(tangents)){
                        # covariance between sample/target points and 
                        # directional derivatives
                        dK <- covariance_matrix_d1(x, tangents, model, dip,
                                                   strike, covariance)
                        ddK <- covariance_matrix_d2(tangents, model, dip,
                                                    strike, covariance)
                        # regularization to avoid singular matrix
                        ddK <- ddK + diag(1e-12, nrow(ddK), ncol(ddK))
                        # trend derivative
                        if(!is.null(trend)) 
                                dT <- trend_matrix_d1(tangents, trend, 
                                                      dip, strike)
                }
                Ntrend <- dim(TR)[2]
                        
                
                # building system
                Kaug <- rbind(
                        cbind(K, dK, TR),
                        cbind(t(dK), ddK, dT),
                        cbind(t(TR), t(dT), matrix(0,Ntrend,Ntrend))
                        )
                # K0aug <- rbind(K0, dK0, t(TR0)) ###
                
                # solving system
                # Kaug <- Matrix(Kaug, dimnames = list(NULL,NULL))
                # K0aug <- Matrix(K0aug) ###
                Kaug <- Matrix(Kaug)
                Kinv <- solve(Kaug)
                
                ## result
                ydata <- getData(y)
                wlist <- list()
                
                # sample values
                xval <- rbind(
                        as.matrix(getData(x)[na_ids,value]),
                        matrix(0, Ntang*2, 1)
                )
                
                # slice loop
                for(i in seq(Nslice)){
                        
                        # slice ID
                        slid <- seq((i - 1) * maxgrid + 1,
                                    min(Ngrid, i * maxgrid))
                        ytemp <- suppressWarnings(y[slid,])
                        
                        # target covariance matrix
                        K0 <- covariance_matrix(x, ytemp, model = model,
                                                covariance = covariance)
                        
                        # trend matrix
                        if(is.null(trend)){
                                TR0 <- matrix(0,length(slid),0)
                        }else{
                                TR0 <- trend_matrix(ytemp, trend)
                        }
                        
                        # tangents
                        dK0 <- matrix(0, 0, length(slid))
                        if(!is.null(tangents)){
                                dK0 <- t(covariance_matrix_d1(ytemp, 
                                                              tangents, model, 
                                                              dip, strike, 
                                                              covariance)) 
                        }
                        
                        # full covariance matrix
                        K0aug <- rbind(K0, dK0, t(TR0))
                        K0aug <- Matrix(K0aug)
                        
                        # solution
                        krigsol <- Kinv %*% K0aug
                
                        # back to normal matrix
                        krigsol <- as.matrix(krigsol) 
                        K0aug <- as.matrix(K0aug)
                
                        # weights
                        lambda <- krigsol[1:(Ndata+2*Ntang),,drop=F]
                
                        # mean
                        if(is.null(mean)){
                                meanval <- matrix(0,length(slid),1)
                        }else{
                                meanval <- as.matrix(ydata[slid,mean])
                        }
                
                        # estimate
                        est <- t(lambda) %*% xval + # probabilistic part
                                matrix(1 - colSums(lambda[1:Ndata,]), 
                                       length(slid), 1) * meanval # mean
                        ydata[slid,to] <- as.numeric(est)
                
                        # kriging variance
                        if(is.null(trend)){
                                ydata[slid, paste0(to,".krigvar")] <-
                                        totvar - colSums(krigsol * K0aug)
                        }else{
                                ydata[slid, paste0(to,".krigvar")] <-
                                        colSums(krigsol * K0aug)
                        }
                
                        # kriging weights
                        if(weights){
                                lnames <- paste0("data", seq(Ndata)[na_ids])
                                if(!is.null(tangents)){
                                        lnames <- c(
                                                lnames,
                                                paste0("dip", seq(Ntang)),
                                                paste0("strike", seq(Ntang))
                                        )
                                }
                                rownames(lambda) <- lnames
                                wlist2 <- apply(lambda, 2, function(rw) list(rw))
                                wlist2 <- lapply(wlist2, unlist)
                                names(wlist2) <- NULL
                                wlist <- c(wlist, wlist2)
                                
                        }
                
                }
                
                ## end
                if(weights) ydata[[paste0(to,".weights")]] <- I(wlist)
                y@data <- ydata
                return(y)
        }
)

setMethod(
        f = "krig3D",
        signature = c("ANY","missing"),
        definition = function(x, y, model, value, CVpart = 5, 
                              nugget = 0, mean = NULL, 
                              trend = NULL, verbose = TRUE, 
                              tangents = NULL, dip = "Dip", strike = "Strike"){
                # Cross validation when y is missing
                require(Matrix)
                
                # setup
                xdata <- getData(x)
                if(!(value %in% colnames(xdata))){
                        stop(paste0("Inexistent value '",value,"'\n"))
                }
                
                # filtering NAs
                na_ids <- !is.na(getData(x)[,value]) 
                
                # counting dimensions
                Ndata <- nrow(x[na_ids,])
                Ngrid <- Ndata
                Ntang <- ifelse(is.null(tangents), 0, nrow(tangents))
                Ntrend <- 0
                
                # parsing nugget
                if(length(nugget) == 1){
                        if(class(nugget) == "character"){
                                nugget <- xdata[,nugget]
                                nugget <- diag(nugget[na_ids])
                        }else{
                                nugget <- diag(nugget,Ndata,Ndata)
                        }
                }else if(is.null(dim(nugget)) || min(dim(nugget)) == 1){
                        nugget <- diag(nugget[na_ids])
                }
                
                # parsing covariance model
                if(length(model) > 1){
                        totvar <- sum(sapply(model, function(m) m@contribution))
                }else{
                        totvar <- model@contribution
                }
                covariance <- is.null(trend)
                
                # parsing CVpart
                if(is.character(CVpart)){
                        CVpart <- xdata[, CVpart]
                }else if(length(CVpart) == 1){
                        CVpart <- sample(CVpart, nrow(x), replace = T)
                }
                
                ## building system
                # covariance matrix
                K <- covariance_matrix(x, x, model = model, covariance = covariance)
                
                # trend matrix
                if(is.null(trend)){
                        TR <- matrix(0,Ndata,0)
                }else{
                        trend <- as.formula(toupper(trend))
                        TR <- model.matrix(trend, getCoords(x, "data.frame"))
                        nugget <- -nugget
                        Ntrend <- dim(TR)[2]
                }
                
                # tangents
                if(is.null(tangents)){
                        dK <- matrix(0,Ndata,0)
                        ddK <- matrix(0,0,0)
                        dT <- matrix(0,0,Ntrend)
                }else{
                        # fazer
                }
                
                # building system
                Kaug <- rbind(
                        cbind(K, dK, TR),
                        cbind(t(dK), ddK, dT),
                        cbind(t(TR), t(dT), matrix(0,Ntrend,Ntrend))
                )
                K0aug <- rbind(K, t(dK), t(TR))
                nugget <- rbind(
                        cbind(nugget, matrix(0,Ndata,Ntrend+2*Ntang)),
                        matrix(0,Ntrend+2*Ntang,Ndata+Ntrend+2*Ntang)
                )
                
                # cross validation
                if(verbose) cat("Performing cross validation ...\n")
                nugget <- Matrix(nugget, dimnames = list(NULL,NULL))
                Kaug <- Matrix(Kaug, dimnames = list(NULL,NULL))
                K0aug <- Matrix(K0aug)
                CVcalc <- rep(NA,Ndata)
                for(item in unique(CVpart)){
                        ids <- which(CVpart == item)
                        krigsol <- solve(Kaug[-ids,-ids] + nugget[-ids,-ids], 
                                         K0aug[-ids,ids])
                        xval <- as.matrix(xdata[na_ids,value,drop=F])
                        xval <- xval[-ids,,drop=F]
                        if(is.null(tangents)){
                                lambda <- krigsol[1:(Ndata-length(ids)),,drop=F]
                                
                        }else{
                                # fazer
                        }
                        # mean
                        if(is.null(mean)){
                                meanval <- matrix(0,length(ids),1)
                        }else{
                                meanval <- as.matrix(xdata[ids,mean])
                        }
                        # Z calc
                        CVcalc[ids] <- t(lambda) %*% xval + # probabilistic part
                                matrix(1 - colSums(lambda), length(ids), 1) * 
                                      meanval # mean
                }
                
                # result
                xval <- as.matrix(xdata[na_ids,value,drop=F])
                xdata[na_ids, paste0(value,".CV")] <- CVcalc
                x@data <- xdata
                return(x)
        }
)

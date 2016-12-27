#### Gaussian Process class ####
GP <- setClass(
        "GP",
        slots = c(data = "spatial3DDataFrame",
                  tangents = "points3DDataFrame",
                  model = "list",
                  mean = "numeric",
                  trend = "character",
                  beta = "matrix",
                  likelihood = "numeric",
                  pre_comp = "list"),
        validity = function(object) {
                if(!all(rapply(object@model,class) == "covarianceStructure3D"))
                        stop("Invalid covariance object")
        }
)

#### initialization ####
setMethod(
        f = "initialize",
        signature = "GP",
        definition = function(.Object, data, model, value, nugget, 
                              mean = NULL, trend = NULL, weights = NULL,
                              tangents = NULL, reg.t = "reg"){
                require(Matrix)
                
                # covariance model
                if(length(model) == 1 & class(model) != "list")
                        .Object@model <- list(model)
                else
                        .Object@model <- model
                
                # value
                if(length(value) == 1 & class(value) == "character")
                        data["value"] <- data[value]
                else
                        data["value"] <- value
                yval <- data[["value"]]
                
                # nugget
                if(length(nugget) == 1){
                        if(class(nugget) == "character"){
                                data["nugget"] <- data[nugget]
                        }else{
                                data["nugget"] <- rep(nugget, nrow(data))
                        }
                }
                else{
                        data["nugget"] <- nugget
                }
                nugget <- data[["nugget"]]
                
                # weights
                if(is.null(weights)) weights <- rep(1, nrow(data))
                data["weights"] <- weights
                
                # regularization for tangents
                if(!is.null(tangents) && nrow(tangents) > 0){
                        if(length(reg.t) == 1){
                                if(class(reg.t) == "character"){
                                        tangents["reg"] <- tangents[reg.t]
                                }else{
                                        tangents["reg"] <- rep(reg.t, nrow(tangents))
                                }
                        }
                        else{
                                tangents["reg"] <- reg.t
                        }
                        reg.t <- tangents[["reg"]]
                }
                
                # mean
                if(is.null(mean)) mean <- mean(yval)
                
                # data
                .Object@data <- data[c("value", "nugget", "weights")]
                .Object@tangents <- as(tangents, "points3DDataFrame")
                
                # trend
                if(is.null(trend) | length(trend) == 0){
                        TR <- matrix(0,nrow(data),0)
                }else{
                        TR <- trend_matrix(data, trend)
                        if(!is.null(tangents) && nrow(tangents) > 0){
                                TR <- rbind(
                                        TR,
                                        trend_matrix_d1(tangents, trend)
                                )
                        }
                        mean <- 0
                }
                Ntrend <- dim(TR)[2]
                
                # covariances
                Ntang <- 0
                K <- covariance_matrix(data, data, model, T) + 
                        diag(nugget / (weights + 1e-6), 
                             length(nugget), length(nugget))
                if(!is.null(tangents) && nrow(tangents) > 0){
                        K1 <- covariance_matrix_d1(data, tangents, model, T)
                        K2 <- covariance_matrix_d2(tangents, model, T)
                        Ntang <- nrow(K2)
                        # regularization
                        K2 <- K2 + diag(reg.t, Ntang, Ntang)
                        # final matrix
                        K <- rbind(
                                cbind(K, K1),
                                cbind(t(K1), K2)
                        )
                }
                
                
                # pre-computations
                .Object@pre_comp <- list()
                L <- t(chol(Matrix(K)))  # makes L lower triangular
                .Object@mean <- mean
                yval <- c(yval - mean, rep(0, Ntang))
                LiY <- solve(L, Matrix(yval, length(yval), 1))
                w_value <- solve(t(L), LiY)
                .Object@pre_comp$w_value <- as.numeric(w_value)
                .Object@pre_comp$w_var <- L
                if(Ntrend > 0){
                        HLi <- solve(L, TR)
                        A <- t(HLi) %*% HLi
                        b1 <- t(TR) %*% w_value
                        beta <- solve(A, b1)
                        .Object@beta <- as.matrix(beta)
                        .Object@trend <- trend
                        .Object@pre_comp$w_trend <- HLi
                }
                
                # likelihood
                dt <- 2*sum(diag(L)^2)
                .Object@likelihood <- -0.5 * dt -
                        0.5 * sum(yval * w_value) -
                        0.5 * length(yval) * log(2*pi)
                if(Ntrend > 0){
                        LiYH <- t(LiY) %*% HLi
                        tmp <- LiYH %*% solve(A, t(LiYH))
                        .Object@likelihood <- .Object@likelihood + 
                                0.5 * as.numeric(tmp) + 
                                Ntrend * log(2*pi) - 
                                0.5 * determinant(A)$modulus
                }
                
                # end
                validObject(.Object)
                return(.Object)
        }
)

#### show ####
setMethod(
        f = "show",
        signature = "GP",
        definition = function(object){
                # display
                cat("Object of class ", class(object), "\n", sep = "")
                cat("Data points:", nrow(object@data), "\n")
                cat("Tangent points:", nrow(object@tangents), "\n")
                if(length(object@trend) == 0)
                        cat("Global mean:", object@mean, "\n")
                else{
                        cat("Trend:", object@trend, "\n")
                        print(object@beta)
                }
                cat("Log-likelihood:", object@likelihood, "\n")
        }
)

#### predict ####
setMethod(
        f = "predict",
        signature = "GP",
        definition = function(object, target, to = "value", output.var = F){
                require(Matrix)
                
                # pre processing
                w_var <- object@pre_comp$w_var
                w_value <- object@pre_comp$w_value
                
                # trend
                if(length(object@trend) > 0){
                        TRtarget <- trend_matrix(target, object@trend)
                        TRdata <- trend_matrix(object@data, object@trend)
                        w_tr <- object@pre_comp$w_trend
                        beta <- object@beta
                }
                
                # slicing target to save memory
                maxgrid <- 1000 # optimize this
                Ngrid <- nrow(target)
                Nslice <- ceiling(Ngrid / maxgrid)
                t2 <- pointify(target)
                
                for(i in seq(Nslice)){
                        
                        # slice ID
                        slid <- seq((i - 1) * maxgrid + 1,
                                    min(Ngrid, i * maxgrid))
                        ttemp <- t2[slid,]
                        
                        # covariances
                        Ntang <- nrow(object@tangents)
                        Ktarget <- covariance_matrix(ttemp, object@data, 
                                                     object@model, T)
                        if(Ntang > 0){
                                K1 <- covariance_matrix_d1(ttemp, 
                                                           object@tangents, 
                                                           object@model, T)
                                Ktarget <- cbind(Ktarget, K1)
                        }
                        
                        # prediction
                        # residuals
                        
                        pred <- apply(Ktarget, 1, function(rw){
                                sum(rw * w_value)
                        }) + object@mean
                        # trend
                        if(length(object@trend) > 0){
                                LinvK <- solve(w_var, t(Ktarget))
                                R <- t(TRtarget[slid,]) - t(w_tr) %*% LinvK
                                pred <- pred + t(R) %*% beta
                        }
                        target[slid, to] <- pred
                        
                        # variance
                        if(output.var){
                                tot_var <- sum(sapply(object@model, 
                                                      function(m) m@contribution))
                                if(length(object@trend) > 0){
                                        pred_var <- colSums(LinvK^2)
                                        pred_var[pred_var > tot_var] <- tot_var
                                        pred_var <- tot_var - pred_var
                                        tr_var <- colSums(
                                                R * (solve(t(w_tr) %*% w_tr, R))
                                                )
                                        pred_var <- pred_var + tr_var
                                }
                                else{
                                        pred_var <- colSums(solve(w_var, t(Ktarget))^2)
                                        pred_var[pred_var > tot_var] <- tot_var
                                        pred_var <- tot_var - pred_var
                                }
                                target[slid, paste0(to, ".var")] <- pred_var
                        }
                
                }
                
                # output
                return(target)
        }
)

#### log-likelihood ####
setMethod(
        f = "logLik",
        signature = "GP",
        definition = function(object){
                return(object@likelihood)
        }
)

#### fit ####
setMethod(
        f = "fit",
        signature = "GP",
        definition = function(object,
                              midrange = F, minrange = F,
                              azimuth = F, dip = F, rake = F,
                              power = F, nugget = F, 
                              nugget.fix = numeric(),
                              seed = runif(1, 1, 10000)){
                require(GA)
                
                # setup
                structures <- sapply(object@model, function(x) x@type)
                Nstruct <- length(structures)
                Ndata <- nrow(object@data)
                data_var <- var(object@data[["value"]])
                data_box <- boundingBox(object@data)
                data_rbase <- sqrt(sum(data_box[1,] - data_box[2,])^2)
                data_nugget <- object@data[["nugget"]]
                
                # optimization limits and starting point
                opt_min <- opt_max <- numeric(Nstruct * 8 + 1)
                xstart <- matrix(0, 1, Nstruct * 8 + 1)
                for(i in 1:Nstruct){
                        # contribution
                        opt_min[(i-1)*8+1] <- data_var/1000
                        opt_max[(i-1)*8+1] <- data_var*2
                        xstart[(i-1)*8+1] <- object@model[[i]]@contribution
                        
                        # maxrange
                        opt_min[(i-1)*8+2] <- data_rbase/1000
                        opt_max[(i-1)*8+2] <- data_rbase
                        xstart[(i-1)*8+2] <- object@model[[i]]@maxrange
                        
                        # midrange (multiple of maxrange)
                        if(midrange)
                                opt_min[(i-1)*8+3] <- 0.01
                        else
                                opt_min[(i-1)*8+3] <- 1
                        opt_max[(i-1)*8+3] <- 1
                        xstart[(i-1)*8+3] <- object@model[[i]]@midrange /
                                object@model[[i]]@maxrange
                        
                        # minrange(multiple of midrange)
                        if(minrange)
                                opt_min[(i-1)*Nstruct+4] <- 0.01
                        else
                                opt_min[(i-1)*Nstruct+4] <- 1
                        opt_max[(i-1)*8+4] <- 1
                        xstart[(i-1)*8+4] <- object@model[[i]]@minrange / 
                                object@model[[i]]@midrange
                        
                        # azimuth
                        opt_min[(i-1)*8+5] <- 0
                        if(azimuth)
                                opt_max[(i-1)*8+5] <- 360
                        else
                                opt_max[(i-1)*8+5] <- 0
                        xstart[(i-1)*8+5] <- object@model[[i]]@azimuth
                        
                        # dip
                        opt_min[(i-1)*8+6] <- 0
                        if(dip)
                                opt_max[(i-1)*8+6] <- 90
                        else
                                opt_max[(i-1)*8+6] <- 0
                        xstart[(i-1)*8+6] <- object@model[[i]]@dip
                        
                        # rake
                        opt_min[(i-1)*8+7] <- 0
                        if(rake)
                                opt_max[(i-1)*8+7] <- 90
                        else
                                opt_max[(i-1)*8+7] <- 0
                        xstart[(i-1)*8+7] <- object@model[[i]]@rake
                        
                        # power
                        if(power){
                                opt_min[(i-1)*8+8] <- 0.1
                                opt_max[(i-1)*8+8] <- 3
                        }
                        else{
                                opt_min[(i-1)*8+8] <- 1
                                opt_max[(i-1)*8+8] <- 1
                        }
                        xstart[(i-1)*8+8] <- object@model[[i]]@power
                }
                
                # nugget
                if(nugget){
                        opt_min[Nstruct * 8 + 1] <- data_var/1000
                        opt_max[Nstruct * 8 + 1] <- data_var*2
                }
                else{
                        opt_min[Nstruct * 8 + 1] <- 0 # not used
                        opt_max[Nstruct * 8 + 1] <- 0 # not used
                }
                xstart[Nstruct * 8 + 1] <- mean(data_nugget)
                
                # conforming starting point to limits
                xstart[xstart < opt_min] <- opt_min[xstart < opt_min]
                xstart[xstart > opt_max] <- opt_max[xstart > opt_max]
                
                # fitness function
                makeGP <- function(x, finished = F){
                        # covariance model
                        m <- vector("list", Nstruct)
                        for(i in 1:Nstruct){
                                m[[i]] <- covarianceStructure3D(
                                        type = structures[i],
                                        contribution = x[(i-1)*8+1],
                                        maxrange = x[(i-1)*8+2],
                                        midrange = x[(i-1)*8+2] * 
                                                x[(i-1)*8+3],
                                        minrange = x[(i-1)*8+2] * 
                                                x[(i-1)*8+3] * 
                                                x[(i-1)*8+4],
                                        azimuth = x[(i-1)*8+5],
                                        dip = x[(i-1)*8+6],
                                        rake = x[(i-1)*8+7],
                                        power = x[(i-1)*8+8]
                                )
                        }
                        # temporary GP
                        if(nugget){ 
                                # fit a constant nugget model
                                tmpnug <- x[Nstruct * 8 + 1] 
                        }
                        else{ 
                                # use values as given
                                tmpnug <- data_nugget
                        }
                        # points with fixed nugget
                        tmpnug[nugget.fix] <- data_nugget[nugget.fix]
                        # GP
                        tmpgp <- GP(
                                data = object@data,
                                model = m,
                                value = "value",
                                nugget = tmpnug,
                                mean = object@mean,
                                trend = object@trend,
                                tangents = object@tangents,
                                weights = object@data[["weights"]]
                        )
                        # output
                        if(finished)
                                return(tmpgp)
                        else
                                return(tmpgp@likelihood)
                }
                
                # optimization
                opt <- ga(
                        type = "real-valued",
                        fitness = function(x) makeGP(x, F),
                        min = opt_min,
                        max = opt_max,
                        pmutation = 0.5,
                        popSize = 20,
                        run = 20,
                        seed = seed,
                        monitor = F,
                        suggestions = xstart
                )
                
                # update
                sol <- opt@solution
                return(makeGP(sol, T))
        }
)
#### Gaussian Process for Logistic Regression ####
GP_LR <- setClass(
        "GP_LR",
        slots = c(GPs = "list",
                  params = "list")
)

#### initialization ####
setMethod(
        f = "initialize",
        signature = "GP_LR",
        definition = function(.Object, data, value1, value2 = value1,
                              model, nugget,
                              tangents = NULL, reg.v = 1e-9, reg.t = 1e-9){
                # setup
                Ndata <- nrow(data)
                xdata <- getData(data)
                categories <- sort(unique(c(xdata[,value1], xdata[,value2])))
                ncat <- length(categories)
                GPs <- vector("list", ncat)

                # model building
                for(i in seq_along(categories)){
                        # balanced indicators
                        ind <- matrix(-1/ncat, nrow(data), 2) # negative
                        ind[xdata[,value1] == categories[i], 1] <- 1 # positive 1
                        ind[xdata[,value2] == categories[i], 2] <- 1 # positive 2
                        ind <- rowMeans(ind) # contacts get (1 - 1/ncat)/2
                        
                        # heteroscedastic nugget
                        nug <- rep(nugget, length(ind))
                        nug[ind == (1 - 1/ncat)/2] <- reg.v # regularization
                        
                        # GP
                        GPs[[i]] <- GP(data, model, 
                                       ind, nug, 
                                       mean = - 1/ncat,
                                       tangents = tangents,
                                       reg.t = reg.t)
                        GPs[[i]]@data["is.contact"] <- ind == (1 - 1/ncat)/2
                }
                names(GPs) <- make.names(categories)
                .Object@GPs <- GPs
                
                # output
                .Object@params$nugget <- nugget
                .Object@params$reg.v <- reg.v
                .Object@params$reg.t <- reg.t
                
                return(.Object)
        }
)

#### show ####
setMethod(
        f = "show",
        signature = "GP_LR",
        definition = function(object){
                # display
                cat("Object of class ", class(object), "\n", sep = "")
                cat("Data points:", nrow(object@GPs[[1]]@data), "\n")
                Ntang <- ifelse(is.null(object@GPs[[1]]@tangents), 0, 
                                nrow(object@GPs[[1]]@tangents))
                cat("Tangent points:", Ntang, "\n")
                cat("Number of classes:", length(object@GPs), "\n")
                cat("Log-likelihood:", logLik(object), "\n")
                lab <- names(object@GPs)
                cat("Class labels:\n")
                for(i in seq_along(lab))
                        cat("  ", lab[i], "\n")
        }
)

#### log-likelihood ####
setMethod(
        f = "logLik",
        signature = "GP_LR",
        definition = function(object){
                ll <- sum(sapply(object@GPs, function(gp) logLik(gp)))
                return(ll)
        }
)

#### predict ####
setMethod(
        f = "predict",
        signature = "GP_LR",
        definition = function(object, target, to = "value", output.unc = F){
                # setup
                ncat <- length(object@GPs)
                indmat <- matrix(0, nrow(target), ncat)
                varmat <- matrix(0, nrow(target), ncat)
                ydata <- getData(target)
                # strength <- max(abs(object@GPs[[1]]@data[["value"]]))
                
                # prediction
                for(i in seq(ncat)){
                        tmp <- predict(object@GPs[[i]], target,
                                       to = "value", output.var = output.unc)
                        indmat[,i] <- tmp[["value"]]
                        if(output.unc)
                                varmat[,i] <- tmp[["value.var"]]
                }
                
                # unknown class
                indmat <- cbind(indmat, -rowSums(indmat)) # zero sum indicators
                
                # indicators
                # skewed potential calculated from probabilities
                indmat_sk <- t(apply(indmat, 1, function(x){
                        xsort <- sort(x, decreasing = T)
                        x - mean(xsort[1:2])
                }))
                inddf <- data.frame(indmat_sk)
                colnames(inddf) <- 
                        paste0(to,"..", c(names(object@GPs), "Unknown"),".ind")
                ydata[,colnames(inddf)] <- inddf
                
                # probabilities and uncertainty
                if(output.unc){
                        probmat <- indmat
                        Nsamples <- 100
                        entropy <- numeric(nrow(probmat))
                        for(i in seq(nrow(probmat))){
                                # sampling
                                sampmat <- matrix(0, Nsamples, ncat)
                                for(j in seq(ncat)){
                                        sampmat[,j] <- rnorm(Nsamples, 
                                                              indmat[i,j],
                                                              sqrt(varmat[i,j]))
                                }
                                sampmat <- cbind(sampmat, -rowSums(sampmat))
                                
                                # probabilities
                                winner <- apply(sampmat, 1, which.max)
                                prob <- sapply(seq(ncat + 1), 
                                               function(x) mean(winner == x))
                                probmat[i,] <- prob
                                
                                # normalized entropy
                                entropy[i] <- 
                                        - sum(log(prob+1e-9) * prob) / log(ncat + 1)
                        }
                        
                        # colnames(probmat) <- 
                        #         paste0(to,"..", c(names(object@GPs), "Unknown"), 
                        #                ".prob")
                        # ydata[,colnames(probmat)] <- as.data.frame(probmat)
                        ydata[,paste0(to,"..entropy")] <- entropy
                }
                
                # # probabilities and classification
                # calcprob <- function(rw){
                #         # preventing NaNs
                #         rw[rw > 20] <- 20
                #         rw[rw < -20] <- -20
                #         # randomly breaking ties
                #         rw <- rw + rnorm(length(rw), sd = 1e-6)
                #         # probabilities
                #         exp(rw) / sum(exp(rw))
                # }
                # probmat <- t(apply(indmat, 1, calcprob))
                # ids <- apply(probmat, 1, which.max)
                # ydata[,to] <- names(object@GPs)[ids]
                # colnames(probmat) <- paste0(to,"..", names(object@GPs), ".prob")
                # ydata[,colnames(probmat)] <- as.data.frame(probmat)
                # 
                # # indicators
                # # skewed potential calculated from probabilities
                # indmat_sk <- t(apply(probmat, 1, function(x){
                #         y <- log(x)
                #         ysort <- sort(y, decreasing = T)
                #         y - mean(ysort[1:2])
                # }))
                # inddf <- data.frame(indmat_sk)
                # colnames(inddf) <- paste0(to,"..", names(object@GPs),".ind")
                # ydata[,colnames(inddf)] <- inddf
                # 
                # # uncertainty (square root of relative variance)
                # if(output.unc){
                #         model <- object@GPs[[1]]@model
                #         if(class(model) != "list") model <- list(model)
                #         totvar <- sum(sapply(model, function(m) m@contribution))
                #         # entropy <- rowSums(- probmat * log(probmat)) / log(ncat)
                #         rel_var <- apply(varmat, 1, min) / totvar
                #         # phi_rel <- apply(indmat, 1, max) / strength
                #         # phi_exp <- min(1, exp(-0.6931*(phi_rel+1)))
                #         u <- sqrt(rel_var)
                #         ydata[, paste0(to,"..uncertainty")] <- u
                # }
                
                # end
                target@data <- ydata
                return(target)
        }
)

#### fit ####
setMethod(
        f = "fit",
        signature = "GP_LR",
        definition = function(object, maxrange = T, midrange = F, minrange = F,
                              azimuth = F, dip = F, rake = F,
                              power = F, 
                              seed = runif(1, 1, 10000)){
                require(GA)
                
                # setup
                structures <- sapply(object@GPs[[1]]@model, 
                                     function(x) x@type)
                Nstruct <- length(structures)
                Ndata <- nrow(object@GPs[[1]]@data)
                # data_var <- 1 #var(object@data[["value"]])
                data_box <- boundingBox(object@GPs[[1]]@data)
                data_rbase <- sqrt(sum(data_box[1,] - data_box[2,])^2)
                # data_nugget <- object@data[["nugget"]]
                GPs <- length(object@GPs)
                
                
                # optimization limits and starting point
                opt_min <- opt_max <- numeric(Nstruct * 8 + 1)
                xstart <- matrix(0, 1, Nstruct * 8 + 1)
                for(i in 1:Nstruct){
                        # contribution
                        opt_min[(i-1)*8+1] <- 0.1
                        opt_max[(i-1)*8+1] <- 5
                        xstart[(i-1)*8+1] <- object@GPs[[1]]@model[[i]]@contribution
                        
                        # maxrange
                        if(maxrange){
                                opt_min[(i-1)*8+2] <- data_rbase/1000
                                opt_max[(i-1)*8+2] <- data_rbase*5
                        }
                        else{
                                opt_min[(i-1)*8+2] <- object@GPs[[1]]@model[[i]]@maxrange
                                opt_max[(i-1)*8+2] <- object@GPs[[1]]@model[[i]]@maxrange
                        }
                        xstart[(i-1)*8+2] <- object@GPs[[1]]@model[[i]]@maxrange
                        
                        # midrange (multiple of maxrange)
                        if(midrange)
                                opt_min[(i-1)*8+3] <- 0.01
                        else
                                opt_min[(i-1)*8+3] <- 1
                        opt_max[(i-1)*8+3] <- 1
                        xstart[(i-1)*8+3] <- object@GPs[[1]]@model[[i]]@midrange /
                                object@GPs[[1]]@model[[i]]@maxrange
                        
                        # minrange(multiple of midrange)
                        if(minrange)
                                opt_min[(i-1)*8+4] <- 0.01
                        else
                                opt_min[(i-1)*8+4] <- 1
                        opt_max[(i-1)*8+4] <- 1
                        xstart[(i-1)*8+4] <- object@GPs[[1]]@model[[i]]@minrange / 
                                object@GPs[[1]]@model[[i]]@midrange
                        
                        # azimuth
                        opt_min[(i-1)*8+5] <- 0
                        if(azimuth)
                                opt_max[(i-1)*8+5] <- 360
                        else
                                opt_max[(i-1)*8+5] <- 0
                        xstart[(i-1)*8+5] <- object@GPs[[1]]@model[[i]]@azimuth
                        
                        # dip
                        opt_min[(i-1)*8+6] <- 0
                        if(dip)
                                opt_max[(i-1)*8+6] <- 90
                        else
                                opt_max[(i-1)*8+6] <- 0
                        xstart[(i-1)*8+6] <- object@GPs[[1]]@model[[i]]@dip
                        
                        # rake
                        opt_min[(i-1)*8+7] <- 0
                        if(rake)
                                opt_max[(i-1)*8+7] <- 90
                        else
                                opt_max[(i-1)*8+7] <- 0
                        xstart[(i-1)*8+7] <- object@GPs[[1]]@model[[i]]@rake
                        
                        # power
                        if(power){
                                opt_min[(i-1)*8+8] <- 0.1
                                opt_max[(i-1)*8+8] <- 3
                        }
                        else{
                                opt_min[(i-1)*8+8] <- 1
                                opt_max[(i-1)*8+8] <- 1
                        }
                        xstart[(i-1)*8+8] <- object@GPs[[1]]@model[[i]]@power
                }
                
                # nugget
                # if(nugget){
                        opt_min[Nstruct * 8 + 1] <- 0.05
                        opt_max[Nstruct * 8 + 1] <- 2
                # }
                # else{
                #         opt_min[Nstruct * 8 + 1] <- 0 # not used
                #         opt_max[Nstruct * 8 + 1] <- 0 # not used
                # }
                xstart[Nstruct * 8 + 1] <- 1 #mean(data_nugget)
                
                # conforming starting point to limits
                xstart[xstart < opt_min] <- opt_min[xstart < opt_min]
                xstart[xstart > opt_max] <- opt_max[xstart > opt_max]
                
                # fitness function
                makeGP <- function(x, finished = F){
                        # weird bug
                        # x[is.nan(x)] <- opt_max[is.nan(x)]
                        # x[is.na(x)] <- opt_max[is.na(x)]
                        
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
                        # if(nugget){ 
                        #         # fit a constant nugget model
                        # tmpnug <- x[Nstruct * 8 + 1] 
                        # }
                        # else{ 
                        #         # use values as given
                        #         tmpnug <- data_nugget
                        # }
                        # # points with fixed nugget
                        # tmpnug[nugget.fix] <- data_nugget[nugget.fix]
                        
                        # GPs
                        tmpgp <- object
                        for(i in 1:GPs){
                                tmpnug <- tmpgp@GPs[[i]]@data[["nugget"]]
                                contacts <- tmpgp@GPs[[i]]@data[["is.contact"]]
                                tmpnug[!contacts] <- x[Nstruct * 8 + 1]
                                # tmpnug <- object@params$nugget
                                # if(finished) # no nugget on contacts
                                        # tmpnug[contacts] <- object@params$reg.v
                                # else # constant nugget during fitting
                                #         tmpnug[contacts] <- x[Nstruct * 8 + 1]
                                tmpgp@GPs[[i]] <- GP(
                                        data = tmpgp@GPs[[i]]@data,
                                        model = m,
                                        value = "value",
                                        nugget = tmpnug,
                                        mean = tmpgp@GPs[[i]]@mean,
                                        trend = tmpgp@GPs[[i]]@trend,
                                        tangents = tmpgp@GPs[[i]]@tangents
                                        )
                        }
                        # output
                        tmpgp@params$nugget <- x[Nstruct * 8 + 1]
                        if(finished)
                                return(tmpgp)
                        else
                                return(logLik(tmpgp))
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
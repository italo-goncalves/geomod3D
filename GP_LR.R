#### Gaussian Process for Logistic Regression ####
GP_LR <- setClass(
        "GP_LR",
        slots = c(GPs = "list")
)

#### initialization ####
setMethod(
        f = "initialize",
        signature = "GP_LR",
        definition = function(.Object, data, value1, value2 = value1,
                              strength = 1, model, nugget,
                              tangents = NULL, reg.v = 1e-9, reg.t = 1e-9){
                # setup
                Ndata <- nrow(data)
                xdata <- getData(data)
                categories <- sort(unique(c(xdata[,value1], xdata[,value2])))
                ncat <- length(categories)
                GPs <- vector("list", ncat)
                # if(length(nugget) < ncat) 
                #         nugget <- rep(nugget, length.out = ncat)
                
                # model building
                for(i in seq_along(categories)){
                        # indicator
                        ind <- matrix(-strength, nrow(data), 2) # negative
                        ind[xdata[,value1] == categories[i], 1] <- strength # positive 1
                        ind[xdata[,value2] == categories[i], 2] <- strength # positive 2
                        ind <- rowMeans(ind) # contacts get 0
                        nug <- rep(nugget, length(ind))
                        nug[ind == 0] <- reg.v # regularization
                        GPs[[i]] <- GP(data, model, 
                                       ind, nug, 
                                       mean = - strength,
                                       tangents = tangents,
                                       reg.t = reg.t)
                }
                names(GPs) <- make.names(categories)
                .Object@GPs <- GPs
                
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
                indmat <- varmat <- matrix(0, nrow(target), ncat)
                ydata <- getData(target)
                strength <- max(abs(object@GPs[[1]]@data[["value"]]))
                
                # prediction
                for(i in seq(ncat)){
                        tmp <- predict(object@GPs[[i]], target,
                                       to = "value", output.var = output.unc)
                        indmat[,i] <- tmp[["value"]]
                        if(output.unc)
                                varmat[,i] <- tmp[["value.var"]]
                }
                
                # probabilities and classification
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
                ydata[,to] <- names(object@GPs)[ids]
                colnames(probmat) <- paste0(to,"..", names(object@GPs), ".prob")
                ydata[,colnames(probmat)] <- as.data.frame(probmat)
                
                # indicators
                # skewed potential calculated from probabilities
                indmat_sk <- t(apply(probmat, 1, function(x){
                        y <- log(x)
                        ysort <- sort(y, decreasing = T)
                        y - mean(ysort[1:2])
                }))
                inddf <- data.frame(indmat_sk)
                colnames(inddf) <- paste0(to,"..", names(object@GPs),".ind")
                ydata[,colnames(inddf)] <- inddf
                
                # uncertainty (square root of relative variance)
                if(output.unc){
                        model <- object@GPs[[1]]@model
                        if(class(model) != "list") model <- list(model)
                        totvar <- sum(sapply(model, function(m) m@contribution))
                        # entropy <- rowSums(- probmat * log(probmat)) / log(ncat)
                        rel_var <- apply(varmat, 1, min) / totvar
                        # phi_rel <- apply(indmat, 1, max) / strength
                        # phi_exp <- min(1, exp(-0.6931*(phi_rel+1)))
                        u <- sqrt(rel_var)
                        ydata[, paste0(to,"..uncertainty")] <- u
                }
                
                # end
                target@data <- ydata
                return(target)
        }
)
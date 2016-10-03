#### Sparse Pseudo-input Gaussian Process class ####
SPGP <- setClass(
        "SPGP",
        slots = c(pseudo_inputs = "matrix",
                  pseudo_nugget = "numeric"),
        contains = "GP"
        # validity = function(object) {
        #         if(!all(rapply(object@model,class) == "covarianceStructure3D"))
        #                 stop("Invalid covariance object")
        # }
)

#### initialization ####
setMethod(
        f = "initialize",
        signature = "SPGP",
        definition = function(.Object, data, model, value, nugget, 
                              mean, pseudo_inputs, 
                              pseudo_nugget = rep(0, nrow(pseudo_inputs))){
                require(Matrix)
                .Object@pseudo_inputs <- pseudo_inputs
                
                # covariance model
                if(length(model) == 1 & class(model) != "list")
                        .Object@model <- list(model)
                else
                        .Object@model <- model
                
                
                # nugget
                if(length(nugget) == 1){
                        if(class(nugget) == "character"){
                                data["nugget"] <- data[nugget]
                        }else{
                                data["nugget"] <- rep(nugget, nrow(data))
                        }
                }
                else{
                        data[, "nugget"] <- nugget
                }
                if(length(pseudo_nugget) == 1)
                        .Object@pseudo_nugget <- rep(pseudo_nugget, nrow(pseudo_inputs))
                else
                        .Object@pseudo_nugget <- pseudo_nugget
                pseudo_nugget <- .Object@pseudo_nugget
                
                # data
                .Object@value <- value
                .Object@data <- data[c(.Object@value, "nugget")]
                
                # covariances
                ps <- points3DDataFrame(
                        pseudo_inputs, 
                        data.frame(.dummy = rep(NA, nrow(pseudo_inputs)))
                )
                K_NM <- Matrix(covariance_matrix(data, ps, model, T))
                K_M <- Matrix(
                        covariance_matrix(ps, ps, model, T) + 
                                diag(pseudo_nugget)
                )
                Q <- K_NM %*% solve(K_M, t(K_NM))

                # pre-computations
                delta <- K_M[1,1] - diag(as.matrix(Q))
                nugget <- unlist(getData(data["nugget"]))
                deltainv <- Matrix(diag(1/(delta + nugget)))
                B <- K_M + t(K_NM) %*% (deltainv %*% K_NM)
                Binv <- solve(Matrix(B))
                K_Minv <- solve(Matrix(K_M))
                yval <- unlist(getData(data[value]))
                .Object@mean <- mean
                yval <- yval - mean
                w_value <- Binv %*% t(K_NM) %*% deltainv %*% 
                        Matrix(yval, length(yval), 1)
                .Object@w_value <- as.numeric(w_value)
                .Object@w_var <- as.matrix(K_Minv - Binv)
                
                # end
                # validObject(.Object)
                return(.Object)
        }
)

#### show ####
setMethod(
        f = "show",
        signature = "SPGP",
        definition = function(object){
                # display
                cat("Object of class ", class(object), "\n", sep = "")
                cat("Data points:", nrow(object@data), "\n")
                cat("Pseudo-inputs:", nrow(object@pseudo_inputs), "\n")
                cat("Value:", object@value, "\n")
        }
)

#### predict ####
setMethod(
        f = "predict",
        signature = "SPGP",
        definition = function(object, target, output.var = F){

                # covariances
                ps <- points3DDataFrame(
                        object@pseudo_inputs,
                        data.frame(.dummy = rep(NA, 
                                                nrow(object@pseudo_inputs)))
                )
                Ktarget <- covariance_matrix(target, ps, object@model, T)
                
                # mean
                w_value <- object@w_value
                pred_mean <- apply(Ktarget, 1, function(rw){
                        sum(rw * w_value)
                })
                target[object@value] <- pred_mean + object@mean
                
                # variance
                if(output.var){
                        w_var <- object@w_var
                        tot_var <- sum(sapply(object@model, 
                                              function(m) m@contribution))
                        pred_var <- rowSums((Ktarget %*% w_var) * Ktarget)
                        target[paste0(object@value, ".var")] <- 
                                tot_var - pred_var
                }

                # output
                return(target)
        }
)
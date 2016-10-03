#### Gaussian Process class ####
GP <- setClass(
        "GP",
        slots = c(data = "spatial3DDataFrame",
                  tangents = "points3DDataFrame",
                  model = "list",
                  w_value = "numeric",
                  w_var = "matrix",
                  mean = "numeric",
                  likelihood = "numeric"),
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
                              mean = 0, tangents = NULL, reg.t = 1e-9){
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
                yval <- unlist(getData(data["value"]))
                
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
                nugget <- unlist(getData(data["nugget"]))
                
                # data
                .Object@data <- data[c("value", "nugget")]
                .Object@tangents <- as(tangents, "points3DDataFrame")
                
                # covariances
                Ntang <- 0
                K <- covariance_matrix(data, data, model, T) + diag(nugget)
                if(!is.null(tangents)){
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
                # Kinv <- solve(Matrix(K))
                
                
                # pre-computations
                L <- t(chol(Matrix(K)))  # makes L lower triangular
                .Object@mean <- mean
                yval <- c(yval - mean, rep(0, Ntang))
                # w_value <- solve(K, Matrix(yval, length(yval), 1))
                w_value <- solve(t(L)) %*% solve(L, Matrix(yval, length(yval), 1))
                .Object@w_value <- as.numeric(w_value)
                # .Object@w_var <- as.matrix(K)
                .Object@w_var <- as.matrix(L)
                
                # likelihood
                # dt <- determinant(K)$modulus
                dt <- sum(diag(L)^2)
                .Object@likelihood <- -0.5 * dt -
                        0.5 * sum(yval * w_value) -
                        0.5 * length(yval) * log(2*pi)
                
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
                Ntang <- ifelse(is.null(object@tangents), 0, 
                                nrow(object@tangents))
                cat("Tangent points:", Ntang, "\n")
        }
)

#### predict ####
setMethod(
        f = "predict",
        signature = "GP",
        definition = function(object, target, to = "value", output.var = F){
                require(Matrix)
                
                # pre processing
                if(output.var) w_var <- solve(Matrix(object@w_var))
                
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
                        
                        # mean
                        w_value <- object@w_value
                        pred_mean <- apply(Ktarget, 1, function(rw){
                                sum(rw * w_value)
                        })
                        target[slid, to] <- pred_mean + object@mean
                        
                        # variance
                        if(output.var){
                                tot_var <- sum(sapply(object@model, 
                                                      function(m) m@contribution))
                                # pred_var <- rowSums((Ktarget %*% w_var) * Ktarget)
                                pred_var <- colSums((w_var %*% t(Ktarget))^2)
                                pred_var[pred_var > tot_var] <- tot_var
                                target[slid, paste0(to, ".var")] <- 
                                        tot_var - pred_var
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
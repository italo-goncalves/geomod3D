#### Gaussian Process class with eigen decomposition ####
# only works with constant nugget
GP_fast <- setClass(
        "GP_fast",
        slots = c(egval = "numeric",
                  egvec = "matrix"),
        contains = "GP",
        validity = function(object) {
                if(!all(rapply(object@model,class) == "covarianceStructure3D"))
                        stop("Invalid covariance object")
        }
)

#### initialization ####
setMethod(
        f = "initialize",
        signature = "GP_fast",
        definition = function(.Object, data, model, value, nugget, 
                              mean = 0, tangents = NULL){
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
                        data["value"] <- nugget
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
                K <- covariance_matrix(data, data, model, T)
                if(!is.null(tangents)){
                        K1 <- covariance_matrix_d1(data, tangents, model, T)
                        K2 <- covariance_matrix_d2(tangents, model, T)
                        Ntang <- nrow(K2)
                        K2 <- K2 + diag(1e-9, Ntang, Ntang) # regularization
                        K <- rbind(
                                cbind(K, K1),
                                cbind(t(K1), K2)
                        )
                }
                
                # eigen decomposition
                eg <- eigen(K)
                .Object@egval <- eg$values
                .Object@egvec <- eg$vectors
                
                # pre-computations
                nugget <- c(nugget, rep(0, Ntang))
                .Object@mean <- mean
                yval <- c(yval - mean, rep(0, Ntang))
                delta <- Diagonal(x = 1/(eg$values + nugget))
                egvec <- Matrix(eg$vectors)
                w_value <- (egvec %*% delta) %*% 
                        crossprod(Matrix(eg$vectors), 
                                  Matrix(yval, length(yval), 1))
                .Object@w_value <- as.numeric(w_value)
                .Object@w_var <- as.matrix(egvec %*% sqrt(delta))
                
                # likelihood
                .Object@likelihood <- -0.5 * sum(eg$values + nugget) -
                        0.5 * sum(yval * w_value) -
                        0.5 * length(yval) * log(2*pi)
                
                # end
                # validObject(.Object)
                return(.Object)
        }
)

#### predict ####
setMethod(
        f = "predict",
        signature = "GP_fast",
        definition = function(object, target, to = "value", output.var = F){
                
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
                                w_var <- object@w_var
                                tot_var <- sum(sapply(object@model, 
                                                      function(m) m@contribution))
                                pred_var <- rowSums((Ktarget %*% w_var)^2)
                                # pred_var <- apply(Ktarget %*% w_var, 1, 
                                #                   function(rw) sum(rw^2))
                                target[slid, paste0(to, ".var")] <- 
                                        tot_var - pred_var
                        }
                        
                }
                
                # output
                return(target)
        }
)


#### update ####
setMethod(
        f = "update",
        signature = "GP_fast",
        definition = function(object, value = NULL, nugget = NULL){
                require(Matrix)
                
                # update nugget
                if(is.null(nugget)){
                        nugget <- unlist(getData(
                                object@data["nugget"]
                        ))
                }
                object@data["nugget"] <- nugget
                
                # update value
                if(is.null(value)){
                        nugget <- unlist(getData(
                                object@data["value"]
                        ))
                }
                object@data["value"] <- value
                
                # pre-computations
                Ntang <- ifelse(is.null(object@tangents), 0, 
                                        nrow(object@tangents))
                nugget <- c(nugget, rep(0, Ntang))
                value <- c(value - object@mean, rep(0, Ntang))
                delta <- Diagonal(x = 1/(object@egval + nugget))
                egvec <- Matrix(object@egvec)
                w_value <- (egvec %*% delta) %*% 
                        crossprod(egvec, 
                                  Matrix(value, length(value), 1))
                object@w_value <- as.numeric(w_value)
                object@w_var <- as.matrix(egvec %*% sqrt(delta))
                
                # likelihood
                object@likelihood <- -0.5 * sum(object@egval + nugget) -
                        0.5 * sum(value * w_value) -
                        0.5 * length(value) * log(2*pi)
                
                return(object)
        }
)
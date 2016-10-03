#### Sparse Pseudo-input Gaussian Process with 
#### Partially Independent Conditional approximation ####
SPGP_PIC <- setClass(
        "SPGP_PIC",
        slots = c(centers = "matrix",
                  w_cluster = "list",
                  w_var_NN = "list",
                  w_var_NB = "list",
                  w_var_BB = "list"),
        contains = "SPGP",
        validity = function(object) {
                if(!all(rapply(object@model,class) == "covarianceStructure3D"))
                        stop("Invalid covariance object")
        }
)

#### initialization ####
setMethod(
        f = "initialize",
        signature = "SPGP_PIC",
        definition = function(.Object, data, model, value, nugget, 
                              mean, pseudo_inputs, cluster, centers,
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
                
                # clusters
                if(class(cluster) == "character")
                        data["cluster"] <- data[cluster]
                else
                        data["cluster"] <- cluster
                .Object@centers <- centers
                
                # ordering of data according to clusters
                cl <- unlist(getData(data["cluster"]))
                ord <- sort(cl, index.return = T)
                data <- data[ord$ix, ]
                
                # data
                .Object@value <- value
                .Object@data <- data[c(.Object@value, "nugget", "cluster")]
                
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
                # mean
                .Object@mean <- mean
                # inverse covariance
                delta <- matrix(0, nrow(Q), ncol(Q))
                deltainv <- delta
                nugget <- unlist(getData(data["nugget"]))
                for(i in seq(max(cluster))){
                        id <- which(getData(data["cluster"]) == i)
                        delta[id, id] <-
                                covariance_matrix(data[id,], data[id,], 
                                                  model, T) - 
                                as.matrix(Q[id, id])
                        deltainv[id, id] <- as.matrix(
                                solve(Matrix(delta[id, id] + diag(nugget[id])))
                        )
                }
                deltainv <- Matrix(deltainv)
                B <- K_M + t(K_NM) %*% (deltainv %*% K_NM)
                Binv <- solve(Matrix(B))
                K_Minv <- solve(K_M)
                w_var <- deltainv - deltainv %*% K_NM %*% Binv %*% 
                        t(K_NM) %*% deltainv
                .Object@w_var <- as.matrix(w_var)
                .Object@w_var_NN <- .Object@w_var_NB <- 
                        .Object@w_var_BB <- vector("list", max(cluster))
                for(i in seq(max(cluster))){
                        id <- which(getData(data["cluster"]) == i)
                        .Object@w_var_NN[[i]] <- as.matrix(
                                K_Minv %*% (
                                        t(K_NM) %*% w_var %*% K_NM + 
                                        t(K_NM)[,id] %*% w_var[id,id] %*% K_NM[id,] - 
                                        2 * t(K_NM) %*% w_var[,id] %*% K_NM[id,]
                                ) %*% K_Minv
                        )
                        .Object@w_var_NB[[i]] <- as.matrix(
                                K_Minv %*% (t(K_NM) %*% w_var[,id] - 
                                            t(K_NM)[,id] %*% w_var[id,id])
                        )
                        .Object@w_var_BB[[i]] <- as.matrix(w_var[id,id])
                }
                # global weights
                yval <- unlist(getData(data[value])) - mean
                p <- apply(w_var, 1, function(rw) sum(rw * yval))
                
                .Object@w_value <- apply(K_Minv %*% t(K_NM), 1, 
                                         function(rw) sum(rw * p))
                # cluster weights
                w_cluster <- vector("list", max(cluster))
                for(i in seq(max(cluster))){
                        id <- which(getData(data["cluster"]) == i)
                        w_cluster[[i]] <-
                                apply(K_Minv %*% t(K_NM)[,id], 1, 
                                      function(rw) sum(rw * p[id]))
                }
                .Object@w_cluster <- w_cluster
                # cluster inverse variances
                
                # end
                # validObject(.Object)
                return(.Object)
        }
)

#### predict ####
setMethod(
        f = "predict",
        signature = "SPGP_PIC",
        definition = function(object, target, output.var = F){
                require(Matrix)
                
                # setup
                w_var <- object@w_var
                yval <- unlist(getData(object@data[object@value]))
                yval <- yval - object@mean
                p <- apply(w_var, 1, function(rw) sum(rw * yval))
                
                w_var_NN <- object@w_var_NN
                w_var_NB <- object@w_var_NB
                w_var_BB <- object@w_var_BB
                
                # assigning target points to clusters
                d <- vectorized_pdist(getPoints(target, "matrix"), 
                                      object@centers)
                Tcl <- apply(d, 1, which.min)
                Ncl <- unlist(getData(object@data["cluster"]))
                
                # covariance matrix
                ps <- points3DDataFrame(
                        object@pseudo_inputs,
                        data.frame(.dummy = rep(NA, 
                                                nrow(object@pseudo_inputs)))
                )
                Ktarget <- Matrix(
                        covariance_matrix(target, ps, object@model, T)
                )

                # mean and variance
                w_cluster <- object@w_cluster
                w_value <- object@w_value
                
                pred_mean <- numeric(nrow(target))
                pred_var <- pred_mean
                target2 <- pointify(target)
                for(i in seq(max(Tcl))){
                        # cluster indices
                        Tid <- which(Tcl == i)
                        Nid <- which(Ncl == i)
                        # points in cluster (full covariance)
                        KB <- Matrix(covariance_matrix(target2[Tid],
                                                       object@data[Nid,],
                                                       object@model, T))
                        # mean
                        w_i <- w_value - w_cluster[[i]]
                        p_i <- p[Nid]
                        tmp <- apply(Ktarget[Tid,], 1, function(rw){
                                sum(rw * w_i)
                        }) + apply(KB, 1, function(rw){
                                sum(rw * p_i)
                        })
                        pred_mean[Tid] <- tmp
                        # variance
                        if(output.var){
                                v1 <- rowSums((Ktarget[Tid,] %*% w_var_NN[[i]]) * 
                                                      Ktarget[Tid,])
                                v2 <- 2 * rowSums((Ktarget[Tid,] %*% w_var_NB[[i]]) * 
                                                      KB)
                                v3 <- rowSums((KB %*% w_var_BB[[i]]) * KB)
                                pred_var[Tid] <- v1 + v2 + v3
                        }
                }

                # output
                target[object@value] <- pred_mean + object@mean
                if(output.var){
                        tot_var <- sum(sapply(object@model, 
                                              function(m) m@contribution))
                        target[paste0(object@value, ".var")] <-
                                tot_var - pred_var
                }
                return(target)
        }
)
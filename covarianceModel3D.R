#### covarianceModel3D class ####
covarianceModel3D <- setClass(
        "covarianceModel3D",
        slots = c(parameters = "data.frame",
                  .eval = "function",
                  .eval_grad = "function",
                  .asym = "list",
                  .ev_list = "list",
                  .evgrad_list = "list"),
        validity = function(object) {
                return(TRUE)
        }
)

#### initialize ####
setMethod(
        f = "initialize",
        signature = "covarianceModel3D",
        definition = function(.Object, ...){
                require(pdist)
                
                # getting arguments
                arglist <- list(...)
                # checking
                if(!all(rapply(arglist,class) == "covarianceStructure3D")){
                        stop("All arguments must be of 
                             class 'covarianceStructure3D'")
                }
                
                # building data.frame representation
                N <- length(arglist)
                slot_names <- slotNames(arglist[[1]])
                df <- data.frame()
                for(i in seq(N)){
                        for(j in slot_names){
                                df[i,j] <- slot(arglist[[i]], j)
                        }
                }
                colnames(df) <- slot_names
                .Object@parameters <- df
                
                # building function
                .Object@.ev_list <- vector("list", N)
                .Object@.evgrad_list <- vector("list", N)
                .Object@.asym <- vector("list", N)
                for(i in seq(N)){
                        .Object@.asym[[i]] <- asymmetry3D(
                                df$maxrange[i], df$midrange[i],
                                df$minrange[i], df$azimuth[i],
                                df$dip[i], df$rake[i])
                        if(df$type[i] == "gaussian"){
                                .Object@.ev_list[[i]] <- vm_gauss 
                                .Object@.evgrad_list[[i]] <- dvm_gauss 
                        }else if(df$type[i] == "exponential"){
                                .Object@.ev_list[[i]] <- vm_exp 
                                .Object@.evgrad_list[[i]] <- dvm_exp
                        }
                }
                .Object@.eval <- function(v, u){
                        val <- matrix(0,dim(v)[1],dim(u)[1])
                        for(i in seq_along(.Object@.ev_list)){
                                val <- val + 
                                        .Object@parameters$contribution[i] * 
                                        .Object@.ev_list[[i]](
                                                v, u,
                                                .Object@.asym[[i]])
                        }
                        return(val)
                }
                .Object@.eval_grad <- function(v, u){
                        val <- matrix(0,dim(v)[1],dim(u)[1])
                        for(i in seq_along(.Object@.evgrad_list)){
                                val <- val + 
                                        .Object@parameters$contribution[i] * 
                                        .Object@.evgrad_list[[i]](
                                                v, u,
                                                .Object@.asym[[i]])
                        }
                        return(val)
                }
                
                return(.Object)
        }
)

#### show ####
setMethod(
        f = "show",
        signature = "covarianceModel3D",
        definition = function(object){
                cat("Object of class ", class(object), "\n", sep = "")
                show(object@parameters)
        }
)
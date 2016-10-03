df_to_lines <- function(collar, assay, survey = NULL, 
                        holeid = "HOLEID", from = "FROM", to = "TO", 
                        X = "X", Y = "Y", Z = "Z"){
        # collar, assay, survey -> data frames
        # holeid, from, to, X, Y, Z -> column names
        
        # require(dplyr)
        
        # setup
        collar_names <- colnames(collar)
        collar_names <- collar_names[!(collar_names %in% c(X,Y,Z))]
        assay_names <- colnames(assay)
        assay_names <- assay_names[!(assay_names %in% c(holeid,from,to))]
        
        # processing data    
        df <- merge(collar, assay, by=holeid)
        lines_list <- list()
        if(is.null(survey)){
                for(i in seq(nrow(df))){
                        lines_list[[i]] <- line3D(
                                df[i, c(X,Y,Z)] - c(0,0,df[i,from]),
                                df[i, c(X,Y,Z)] - c(0,0,df[i,to]) 
                        )
                }
        }else{
                # fazer survey aqui
                stop("Survey data not supported yet")
        }
        df <- df[, c(collar_names, assay_names)]
        
        # output
        return(lines3DDataFrame(lines_list,df))
}

get_Google_Elevation <- function(lat, long, key){
        # packages
        require(httr)
        require(jsonlite)
        require(dplyr)
        
        # setup
        Ndata <- length(lat)
        base_int <- 40 # max roughly 40 data points per request
        Ncells <- ceiling(Ndata/base_int) 
        baseURL <- "https://maps.googleapis.com/maps/api/elevation/json?locations="
        
        # conversion to character vector
        lat_long_str <- character(Ncells)
        for (i in 1:(Ncells-1)){
                lat_long_str[i] <- 
                        paste(lat[((i-1)*base_int+1):(i*base_int)],
                              long[((i-1)*base_int+1):(i*base_int)],
                              sep = ",", collapse = "|")
        }
        lat_long_str[Ncells] <- paste(lat[-(1:(i*base_int))],
                                      long[-(1:(i*base_int))],
                                      sep = ",", collapse = "|")
        
        # request to Google Elevation API
        elev <- data.frame()
        for (i in 1:Ncells){
                # Sys.sleep(0.1)  # pause
                url <- paste0(baseURL,lat_long_str[i],"&key=",key)
                json <- fromJSON(
                        content(GET(url), as = "text"), 
                        simplifyDataFrame = T)
                cat("Request",i,"of",Ncells,"- Status =",json$status)
                if(json$status == "OK"){
                        gdata <- as.data.frame(json$results)
                        gdata <- tbl_df(data.frame(gdata$location, 
                                                   gdata$elevation, 
                                                   gdata$resolution))
                        elev <- bind_rows(elev, gdata)
                }
        }
        
        # end
        colnames(elev) <- c("latitude","longitude","elevation","resolution")
        return(elev)
}


rescale_coords <- function(coords, 
                           old_range = apply(coords,2,range), 
                           new_range = matrix(rep(c(0,1),3),2,3)){
        # setup
        require(scales)
        coords <- as.matrix(coords)
        
        # rescaling
        new_coords <- sapply(
                1:3,  
                function(x, coords, old_range, new_range) 
                        scales::rescale(coords[,x], 
                                from = old_range[,x], to = new_range[,x]),
                coords = coords, old_range = old_range, new_range = new_range
                )
        return(new_coords)
}

# linedist <- function(P1, P0, Q1, Q0, as = "dist"){
#         
#         P10 <- P1 - P0
#         P10t <- t(P10)
#         Q10 <- Q1 - Q0
#         Q10t <- t(Q10)
#         PQ0 <- P0 - Q0
#         PQ0t <- t(PQ0)
#         
#         # a <- (P1 - P0) %*% t(P1 - P0)
#         # b <- (P1 - P0) %*% t(Q1 - Q0)
#         # c <- (Q1 - Q0) %*% t(Q1 - Q0)
#         # d <- (P1 - P0) %*% t(P0 - Q0)
#         # e <- (Q1 - Q0) %*% t(P0 - Q0)
#         # f <- (P0 - Q0) %*% t(P0 - Q0)
#         
#         a <- P10 %*% P10t
#         b <- P10 %*% Q10t
#         c <- Q10 %*% Q10t
#         d <- P10 %*% PQ0t
#         e <- Q10 %*% PQ0t
#         
#         if(a*c - b^2 == 0){ # two points
#                 if(as == "dist"){
#                         return(sqrt((P1 - Q1) %*% t((P1 - Q1))))
#                 }else if(as == "coords"){
#                         return((rbind(P1,Q1)))
#                 }
#         }else{
#                 M <- matrix(c(a,-b,-b,c),2,2)
#                 K <- matrix(c(d,-e),2,1)
#                 
#                 sc <- max(abs(cbind(M,K)))
#                 M <- M/sc
#                 K <- K/sc
#                 
#                 p <- -solve(M + diag(1e-9,2,2), K)
#                 alpha <- 1e-1
#                 maxit <- 1000
#                 n_old <- 1e15
#                 zerog1 <- F
#                 zerog2 <- F
#                 for(i in seq(maxit)){
#                         if(max(p) >= 1 | min(p) <= 0){
#                                 if(p[1] >= 1 | p[1] <= 0){
#                                         p[1] <- max(0,min(1,p[1]))
#                                         zerog1 <- T
#                                 }
#                                 if(p[2] >= 1 | p[2] <= 0){
#                                         p[2] <- max(0,min(1,p[2]))
#                                         zerog2 <- T
#                                 }
#                                 grad <- 2 * M %*% p + 2 * K
#                                 if(zerog1) grad[1] <- 0
#                                 if(zerog2) grad[2] <- 0
#                                 # if(sum(abs(grad)) == 0) break
#                                 p <- p - (alpha*(maxit-i)/maxit) * grad
#                         }else{
#                                 grad <- 2 * M %*% p + 2 * K
#                                 p <- p - (alpha*(maxit-i)/maxit) * grad
#                         }
#                         n <- sqrt(sum(grad^2))
#                         if(n < 1e-6 | abs(n - n_old) < 1e-6) break
#                         n_old <- n
#                         # cat(p[1], p[2], n,"\n")
#                 }
#                 
#                 point1 <- (1-p[1])*P0 + p[1]*P1
#                 point2 <- (1-p[2])*Q0 + p[2]*Q1
#                 
#                 if(as == "dist"){
#                         return(sqrt(sum(point1-point2)^2))
#                 }else if(as == "coords"){
#                         return((rbind(point1,point2)))
#                 }
#         }
# }
# 
# linedist2 <- function(P1, P0, Q1, Q0, as = "dist"){
#         require(pdist)
#         
#         
#         break_lines <- function(P1,P0,Q1,Q0, parts){
#                 r <- matrix(rep(seq(0,1,1/parts),times=3),parts+1,3)
#                 if(all(P1 == P0)){
#                         Pmat <- matrix(P1,1,3)
#                 }else{
#                         Pmat <- matrix(P0,parts+1,3, byrow = T) + 
#                                 matrix(P1-P0,parts+1,3, byrow = T) * r
#                 }
#                 if(all(Q1 == Q0)){
#                         Qmat <- matrix(Q1,1,3)
#                 }else{
#                         Qmat <- matrix(Q0,parts+1,3, byrow = T) + 
#                                 matrix(Q1-Q0,parts+1,3, byrow = T) * r
#                 }
#                 return(list(Pmat,Qmat))
#         }
#         
#         # 1st pass
#         br <- break_lines(P1,P0,Q1,Q0,100)
#         suppressWarnings(d <- as.matrix(pdist::pdist(br[[1]],br[[2]])))
#         ids <- which(d == min(d), arr.ind = T)
#         Pmat <- br[[1]]
#         Qmat <- br[[2]]
#         
#         # # 2nd pass
#         # P1 <- Pmat[min(1,ids[1,1]+1),,drop=F]
#         # P0 <- Pmat[max(0,ids[1,1]-1),,drop=F]
#         # Q1 <- Qmat[min(1,ids[1,2]+1),,drop=F]
#         # Q0 <- Qmat[max(0,ids[1,2]-1),,drop=F]
#         # 
#         # br <- break_lines(P1,P0,Q1,Q0,3)
#         # suppressWarnings(d <- as.matrix(pdist::pdist(br[[1]],br[[2]])))
#         # ids <- which(d == min(d), arr.ind = T)
#         # Pmat <- br[[1]]
#         # Qmat <- br[[2]]
#         # 
#         # # 3rd pass
#         # P1 <- Pmat[min(1,ids[1,1]+1),,drop=F]
#         # P0 <- Pmat[max(0,ids[1,1]-1),,drop=F]
#         # Q1 <- Qmat[min(1,ids[1,2]+1),,drop=F]
#         # Q0 <- Qmat[max(0,ids[1,2]-1),,drop=F]
#         # 
#         # br <- break_lines(P1,P0,Q1,Q0,10)
#         # suppressWarnings(d <- as.matrix(pdist::pdist(br[[1]],br[[2]])))
#         # ids <- which(d == min(d), arr.ind = T)
#         # Pmat <- br[[1]]
#         # Qmat <- br[[2]]
#         
#         # end
#         if(as == "dist"){
#                 return(min(d))
#         }else if(as == "coords"){
#                 ids <- which(d == min(d), arr.ind = T)
#                 coords <- rbind(Pmat[ids[1,1],,drop=F],
#                                 Qmat[ids[1,2],,drop=F])
#                 return(coords)
#         }
#         
# }

deriv_formula <- function(f, x){
        require(Deriv)
        # get terms
        f2 <- attr(terms(f), "term.labels")
        if(length(f2) == 0) return(as.formula(paste0("~ I(0*",x,") - 1")))
        intercept <- attr(terms(f), "intercept") == 1
        # symbolic to arithmetic
        f2 <- gsub(":", "*", f2)
        f2 <- gsub(":", "*", f2)
        Is <- grep("I(.+)", f2)
        if(length(Is) > 0){
                f2[Is] <- gsub("[I(]", "", f2[Is])
                f2[Is] <- gsub(")$", "", f2[Is])
        }
        # derivative
        f2 <- sapply(f2, function(a) Deriv(a, x))
        # arithmetic to symbolic
        zeros <- grep("^0$", f2)
        if(length(zeros) > 0){
                f2[zeros] <- paste0("0^",seq_along(zeros),"*",x)
        }
        f2 <- gsub("^1$", paste0(x,"^0"), f2)
        f2 <- paste0("I(", f2, ")", collapse = "+")
        if(intercept) f2 <- paste0("I(0*",x,") + ", f2)
        f2 <- paste("~", f2, "-1")
        f2 <- as.formula(f2)
}

# code from http://blog.felixriedel.com/2013/05/pairwise-distances-in-r/
vectorized_pdist <- function(A,B){
        an <- apply(A, 1, function(rvec) crossprod(rvec,rvec))
        bn <- apply(B, 1, function(rvec) crossprod(rvec,rvec))
        
        m <- nrow(A)
        n <- nrow(B)
        
        tmp <- matrix(rep(an, n), nrow=m)
        tmp <- tmp +  matrix(rep(bn, m), nrow=m, byrow=TRUE)
        d2 <- tmp - 2 * tcrossprod(A,B)
        d2[d2 < 0] <- 0
        sqrt(d2)
}

find_color_cont <- function(values, rng = range(values), col){
        require(scales)
        valsc <- rescale(values, from = rng)
        cr <- colour_ramp(col)
        return(cr(valsc))
}

blockchol <- function(A11, A21, A22){
        require(Matrix)
        L11 <- t(chol(A11))
        L21 <- A21 %*% solve(t(L11))
        tmp <- A22 - L21 %*% t(L21)
        L22 <- t(chol(tmp))
        rbind(cbind(L11, Matrix(0, nrow(L11), ncol(L22))), cbind(L21, L22))
}
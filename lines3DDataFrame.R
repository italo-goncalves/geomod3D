#### lines3DDataFrame class ####
lines3DDataFrame <- setClass(
        "lines3DDataFrame",
        contains = "spatial3DDataFrame",
        validity = function(object) {
                if (length(object@coords) !=
                    nrow(object@data)) {
                        stop(
                                paste(
                                        "Number of line3D objects does not",
                                        "match number of observations"
                                )
                        )
                }
                if (!all(rapply(object@coords,class) == "line3D")) {
                        stop(
                                paste(
                                        "Invalid object in",
                                        "lines list. All",
                                        "objects must be of",
                                        "class 'line3D'"
                                )
                        )
                }
        }
)

#### initialization ####
setMethod(
        f = "initialize",
        signature = "lines3DDataFrame",
        definition = function(.Object, lines_list, df){
                .Object@coords <- lines_list
                .Object@data <- df
                # bounding box
                lines_df <- getCoords(.Object,"data.frame")
                names(lines_df) <- rep(c("X","Y","Z"),2)
                lines_df <- rbind(lines_df[,1:3],lines_df[,4:6])
                bbox <- as.matrix(rbind(
                        apply(lines_df,2,min),
                        apply(lines_df,2,max)))
                rownames(bbox) <- c("min","max")
                .Object@bbox <- bbox
                validObject(.Object)
                return(.Object)
        }
)

#### getCoords ####
setMethod(
        f = "getCoords",
        signature = "lines3DDataFrame",
        definition = function(object, as = "list"){
                if(as == "list"){
                        return(object@coords)
                }else if(as == "data.frame"){
                        lines_list <- object@coords
                        df <- matrix(nrow = length(lines_list), ncol = 6)
                        for(i in seq(length(lines_list))){
                                coords <- getCoords(lines_list[[i]])
                                df[i,] <- cbind(coords[1,], coords[2,])
                        }
                        df <- as.data.frame(df)
                        colnames(df) <- c("X.1","Y.1","Z.1","X.2","Y.2","Z.2")
                        return(df)
                }else if(as == "matrix"){
                        lines_list <- object@coords
                        df <- matrix(nrow = length(lines_list), ncol = 6)
                        for(i in seq(length(lines_list))){
                                coords <- getCoords(lines_list[[i]])
                                df[i,] <- cbind(coords[1,], coords[2,])
                        }
                        colnames(df) <- c("X.1","Y.1","Z.1","X.2","Y.2","Z.2")
                        return(df)
                }else{
                        stop("Invalid value for 'as' argument")
                }
        }
)


#### getLength ####
setMethod(
        f = "getLength",
        signature = "lines3DDataFrame",
        definition = function(object){
                lines_list <- getCoords(object)
                return(sapply(lines_list, getLength))
        }
)


#### show ####
setMethod(
        f = "show",
        signature = "lines3DDataFrame",
        definition = function(object){
                # setup
                l <- min(10, nrow(object))
                l3df <- object[seq(l), ]
                coords <- getCoords(l3df,"data.frame")
                df <- getData(l3df)
                # display
                cat("Object of class ", class(object), "\n", sep = "")
                cat(nrow(object), " line segments and ",
                    ncol(object), " attributes\n\n", sep = "")
                cat("Bounding box:\n")
                show(boundingBox(object))
                cat("\nLine segments:\n")
                show(coords)
                cat("\nAttributes:\n")
                show(df)
        }
)

#### pointify ####
setMethod(
        f = "pointify",
        signature = "lines3DDataFrame",
        definition = function(object, locations = c(0.05,0.5,0.95),
                              distance = F){
                if(min(locations) < 0 || max(locations) > 1){
                        stop("locations must contain values between 0 and 1,
                             inclusive")
                }
                # conversion to points3DDataFrame
                r <- length(locations)
                n <- nrow(object)
                d <- numeric(n * r)
                points_list <- vector("list", n * r)
                lines_list <- getCoords(object)
                df <- getData(object)
                for(i in seq_len(n)){
                        points_list[(r*(i-1)+1):(i*r)] <- 
                                pointify(lines_list[[i]], locations)
                        d[(r*(i-1)+1):(i*r)] <- apply(rbind(
                                getLength(lines_list[[i]]) * locations,
                                getLength(lines_list[[i]]) * (1-locations)
                        ), 2, min)
                }
                df <- df[rep(seq(n), each = r),,drop=F]
                if(distance) df[,".dist"] <- d
                return(points3DDataFrame(points_list,df))
        }
)

#### reScale ####
setMethod(
        f = "reScale",
        signature = "lines3DDataFrame",
        definition = function(object, 
                              old_range = boundingBox(object), 
                              new_range = matrix(rep(c(0,1),3),2,3)){
                lines_list <- getCoords(object)
                df <- getData(object)
                lines_list <- lapply(
                        lines_list,
                        function(x) reScale(x, old_range, new_range)
                )
                return(lines3DDataFrame(lines_list,df))
        }
)

#### drawDrillholes ####
setMethod(
        f = "drawDrillholes",
        signature = "lines3DDataFrame",
        definition = function(object, by, values, col, size){
                # pacakges
                require(rgl)
                require(grDevices)
                require(plotrix)
                
                # setup
                object <- mergeSegments(object, by)
                df <- getData(object)
                lines_list <- getCoords(object)
                if(!any(colnames(df) %in% by)){
                        stop("Invalid attribute")
                }
                N <- nrow(object)
                if(length(size) < N) size <- rep(size, length.out = N)
                objval <- df[,by]
                
                # pallette
                if(class(objval) == "numeric"){ # continuous variable
                        colorsc <- find_color_cont(objval, 
                                                   rng = range(values),
                                                   col = col)
                }else{ # categorical variable
                        names(col) <- values
                        colorsc <- col[objval]
                }
                # plotting
                for(i in seq(N)){
                        shade3d(
                                cylinder3d(getCoords(lines_list[[i]]),
                                           sides = 16, radius = size[i]/2, 
                                           closed = -2),
                                col = colorsc[i])
                }
        }
)

#### drawHoleID ####
setMethod(
        f = "drawHoleID",
        signature = "lines3DDataFrame",
        definition = function(object, cex = 1){
                # setup
                coords <- getCoords(object, "data.frame")
                df <- cbind(coords, getData(object))
                holeids <- unique(df$HOLEID)
                N <- length(holeids)
                coords2 <- as.data.frame(matrix(NA, N, 3))
                
                # highest Z
                for(i in seq(N)){
                        tmp <- df[df$HOLEID == holeids[i], ]
                        pos <- which(tmp$Z.1 == max(tmp$Z.1))
                        coords2[i, ] <- tmp[pos, 1:3]
                }
                
                # plotting
                text3d(coords2, texts = holeids, cex = cex, adj = c(0.5, 0))
                
        }
)

#### mergeSegments ####
setMethod(
        f = "mergeSegments",
        signature = "lines3DDataFrame",
        definition = function(object, by, 
                              keep = colnames(getData(object))){
                require(tidyr)
                require(dplyr)
                
                # setup
                # if(!(categories %in% c("proportional"))){
                #         stop("Invalid value for 'categories' parameter")
                # }
                line_lengths <- getLength(object)
                line_coords <- getCoords(object, "matrix")
                lines_list <- getCoords(object)
                df <- getData(object)
                if(by == "HOLEID"){
                        by <- ".HOLEID"
                        df[,by] <- df[,"HOLEID"]
                }
                if(any(keep %in% c("HOLEID", by))){
                        keep <- keep[-which(keep %in% c("HOLEID", by))]
                }
                df <- df[,c("HOLEID", by, keep)]
                Nlines <- nrow(line_coords)
                directions <- (line_coords[,1:3] - line_coords[,4:6]) / 
                        line_lengths
                
                # finding mergeable segments
                # condition 1 - segments sharing a point
                coord_from <- line_coords[2:Nlines,1:3]
                coord_to <- line_coords[1:(Nlines-1),4:6]
                start_end <- apply(coord_from - coord_to, 1, 
                                   function(x) all(x==0))
                # condition 2 - parallelism
                dir_from <- directions[2:Nlines,]
                dir_to <- directions[1:(Nlines-1),]
                parallel <- apply(dir_from - dir_to, 1,
                                  function(x) all(x==0))
                # condition 3 - same value in "by" column
                val_to <- df[1:(Nlines-1),by]
                val_from <- df[2:Nlines,by]
                same_value <- val_from == val_to
                # final vector
                mergeable <- start_end & parallel & same_value
                
                # building merged object
                # find contiguous mergeable segments
                merge_ids <- split(seq(Nlines), diffinv(!mergeable))
                # number of rows in new object
                Nlines2 <- length(merge_ids)
                # new data frame
                df2 <- df[,-which(colnames(df) %in% c("HOLEID",by)),drop=F]
                cols_to_pass <- colnames(df2)
                df2[,".line"] <- seq(nrow(df2)) # unique ID to avoid spread error
                for(i in cols_to_pass){
                        if(class(df2[,i]) != "numeric"){
                                df2[,i] <- as.character(df2[,i])
                                # special characters in column names may give
                                # error here
                                df2[,i] <- paste0(i,"..", 
                                                  make.names(df2[,i]))
                                df2[,".val"] <- 1
                                df2 <- tidyr::spread_(data = df2,
                                                      key_col = i,
                                                      value_col = ".val",
                                                      fill = 0)
                                df2 <- dplyr::arrange_(df2, ".line")
                        }        
                }
                df2 <- df2[,-which(colnames(df2)==".line")]
                df2 <- cbind(df[,c("HOLEID",by),drop=F], df2)
                # averaging values
                lines_list2 <- vector("list", Nlines2)
                df3 <- data.frame(matrix(NA,Nlines2,ncol(df2)))
                colnames(df3) <- colnames(df2)
                for(i in seq_along(merge_ids)){
                        l <- length(merge_ids[[i]])
                        id_start <- merge_ids[[i]][1]
                        id_end <- merge_ids[[i]][l]
                        lines_list2[[i]] <- line3D(
                                line_coords[id_start,1:3],
                                line_coords[id_end,4:6]
                        )
                        df3[i,c("HOLEID",by)] <- 
                                df2[id_start,c("HOLEID",by)]
                        if(ncol(df3) > 2){
                                if(l > 1){ # average
                                        w <- line_lengths[id_start:id_end] / 
                                                sum(line_lengths[id_start:id_end])
                                        lweights <- matrix(w, l, ncol(df2)-2)
                                        df3[i,seq(3,ncol(df3))] <-
                                                colSums(
                                                        lweights * 
                                                        df2[id_start:id_end, 
                                                            seq(3,ncol(df3)), 
                                                            drop=F]
                                                )
                                }else{ # just copy (faster)
                                        df3[i,seq(3,ncol(df3))] <-
                                                df2[id_start, seq(3,ncol(df3))]
                                }
                        }
                }
                
                # end
                if(by == ".HOLEID") df3 <- 
                        df3[,-which(colnames(df3) == by),drop=F]
                return(lines3DDataFrame(lines_list2, df3))
        }
)

#### smoothingPenalty ####
# setMethod(
#         f = "smoothingPenalty",
#         signature = "lines3DDataFrame",
#         definition = function(object, model, to){
#                 Nrows <- nrow(object)
#                 penalty <- numeric(Nrows)
#                 for(i in seq(Nrows)){
#                         objp <- pointify(object[i,], seq(0,1,0.01))
#                         pts <- getPoints(objp, "matrix")
#                         g <- variogram3D(pts, model = model)
#                         penalty[i] <- mean(g)
#                 }
#                 object[to] <- penalty
#                 return(object)
#         }
# )

#### getContacts ####
setMethod(
        f = "getContacts",
        signature = "lines3DDataFrame",
        definition = function(object, by){
                # setup
                x <- mergeSegments(object, by)
                x <- pointify(x, c(0,1))
                pts <- getCoords(x,"matrix")
                xdata <- getData(x)
                # finding duplicate indices
                dup <- which(duplicated(pts))
                dup2 <- dup-1
                # building new object
                new_points <- getCoords(x)[dup]
                new_df <- data.frame(xdata[dup2,"HOLEID"], 
                                     xdata[dup2,by], xdata[dup,by])
                colnames(new_df) <- c("HOLEID", paste0(by,c(".up",".down")))
                return(points3DDataFrame(new_points,new_df))
        }
)

#### classify ####
setMethod(
        f = "classify",
        signature = c("lines3DDataFrame", "points3DDataFrame"),
        definition = function(x, y, by, model, nugget, trend = "~ x + y + z", 
                              prob = F, uncertainty = F, indicators = T,
                              verbose = T,
                              ortho = NULL, dip = "Dip", strike = "Strike"){
                # setup
                xdata <- getData(x)
                ydata <- getData(y)
                categories <- unique(xdata[,by])
                ncat <- length(categories)
                indmat <- matrix(NA, nrow(y), ncat) # indicators
                varmat <- matrix(NA, nrow(y), ncat) # kriging variance
                # indicator kriging
                for(item in categories){
                        x[".ind"] <- -1 # negative
                        x[getData(x[by]) == item, ".ind"] <- 1 # positive
                        xmerge <- mergeSegments(x, by = ".ind", 
                                                keep = character(0))
                        xmerge[".nugget"] <- 1
                        contacts <- getContacts(xmerge, by = ".ind")
                        contacts[".nugget"] <- 0 # no regularization on contacts
                        contacts[".ind"] <- 0
                        xp <- pointify(xmerge, 
                                       c(0.01,0.05,seq(0.1,0.9,0.2),0.95,0.99))
                        pointdata <- bindPoints(xp, contacts)
                        pnugget <- getData(pointdata[".nugget"]) * nugget
                        pnugget <- unlist(pnugget)
                        if(verbose){
                                cat("Classifying ", by, ": " ,item, " ... ", 
                                    sep = "")
                        }
                        tempgr <- krig3D(pointdata, y, model = model,
                                         value = ".ind", nugget = pnugget,
                                         trend = trend, l1 = l1,
                                         krigvar = T, verbose = verbose,
                                         ortho = ortho, dip = dip, 
                                         strike = strike)
                        indmat[, which(categories == item)] <-
                                unlist(getData(tempgr[".ind"]))
                        varmat[, which(categories == item)] <-
                                unlist(getData(tempgr[".ind.krigvar"]))
                }
                # classification (one vs all)
                code_matrix <- matrix(-1,ncat,ncat) + diag(2,ncat,ncat)
                indmat <- indmat %*% code_matrix
                probmat <- apply(indmat, 1, function(rw){
                        # randomly breaking ties
                        rw <- rw + rnorm(length(rw), sd = 1e-6) 
                        exp(rw) / sum(exp(rw))
                })
                probmat <- t(probmat)
                ids <- apply(probmat, 1, which.max)
                # output
                ydata[,by] <- categories[ids]
                if(prob){
                        colnames(probmat) <- paste0(by,"..",
                                                    make.names(categories),
                                                    ".prob")
                        ydata[,colnames(probmat)] <- as.data.frame(probmat)
                }
                if(indicators){
                        # distance from border calculated from probabilities
                        indmat <- t(apply(probmat, 1, function(x){
                                y <- log(x)
                                ysort <- sort(y, decreasing = T)
                                y - mean(ysort[1:2])
                        }))
                        inddf <- data.frame(indmat)
                        colnames(inddf) <- paste0(by,"..",
                                                  make.names(categories),".ind")
                        ydata[,colnames(inddf)] <- inddf
                }
                if(uncertainty){
                        entropy <- -rowSums(probmat * log(probmat))/log(ncat)
                        avgvar <- rowMeans(varmat)
                        u <- sqrt((1+entropy)*(1+avgvar)) -1
                        ydata[, paste0(by,"..uncertainty")] <- u
                }
                y@data <- ydata
                return(y)
        }
)

setMethod(
        f = "classify",
        signature = c("lines3DDataFrame","missing"),
        definition = function(x, by, model, nugget, trend = "~ x + y + z",  
                              discretization = seq(0.2,0.8,0.2),
                              verbose = TRUE, 
                              ortho = NULL, dip = "Dip", strike = "Strike"){
                # setup
                xdata <- getData(x)
                y <- pointify(mergeSegments(x, by, keep=character(0)), 
                              discretization)
                ydata <- getData(y)
                categories <- unique(xdata[,by])
                ncat <- length(categories)
                # cross-validation
                catcalc <- character(nrow(y))
                for(hole in unique(xdata[,"HOLEID"])){
                        if(verbose) cat("Cross validating hole",hole,"\n")
                        xids <- xdata[,"HOLEID"] != hole
                        yids <- ydata[,"HOLEID"] == hole
                        xtemp <- x[xids,]
                        ytemp <- y[yids,]
                        ytemp <- classify(xtemp,ytemp, by = by, model = model,
                                          nugget = nugget, trend = trend, l1 = l1,
                                          verbose = verbose, ortho = ortho,
                                          dip = dip, strike = strike)
                        catcalc[yids] <- unlist(getData(ytemp[by]))
                }
                cv <- mean(catcalc == ydata[,by])
                return(cv)
        }
)
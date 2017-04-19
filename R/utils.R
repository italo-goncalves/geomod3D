.deriv_formula <- function(f, x){
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
  f2 <- sapply(f2, function(a) Deriv::Deriv(a, x))
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

.find_color_cont <- function(values, rng = range(values), col,
                             na.color = NA){
  valsc <- rescale(values, from = rng)
  cr <- colour_ramp(col, na.color = na.color)
  return(cr(valsc))
}


.onUnload <- function (libpath) {
  library.dynam.unload("geomod3D", libpath)
}

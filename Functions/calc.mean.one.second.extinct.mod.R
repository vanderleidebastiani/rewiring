# Code extracted from the "second.extinct" function to calcultate mean of multiple replicates of extinction sequence and set specific attributes to the object returned
# This way robustness function can be applied to a list of results returned by the function "one.second.extinct.mod"
calc.mean.one.second.extinct.mod <- function(o){
  lengths <- sapply(o, nrow)
  z <- o[[which.max(lengths)]]
  z[, 2:3] <- 0
  for (k in 1:length(o)) {
    nr <- nrow(o[[k]])
    z[1:nr, ] <- z[1:nr, ] + o[[k]]
    rm(nr)
  }
  out <- z/length(o)
  out[, 1] <- 1:max(lengths)
  participant <- attr(o,"exterminated")
  class(out) <- "bipartite"
  attr(out, "exterminated") <- attr(o[[1]],"exterminated")
  return(out)
}

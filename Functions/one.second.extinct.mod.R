# Small change in internal function "one.second.extinct" of bipartite package. The change allows rewiring of interations after each step of extinction
# The arguments are the same of "one.second.extinct" function and additional argumets are described below.
# New arguments:
# rewiring - Logical argument to specify if allow rewiring (default rewiring = FALSE).
# probabilities.rewiring1 - A matrix with probabilities of rewiring, must be the same dimensions of web. See section Methods for details. This matrix is required in step ii of framework  (default probabilities.rewiring1 = NULL).
# probabilities.rewiring2 - A matrix with probabilities of rewiring, must be the same dimensions of web. See section Methods for details. This matrix is required in step iii of framework (default probabilities.rewiring2 = NULL).
# method.rewiring = Type of method used to trial rewiring, partial match to "one.try.single.partner", "multiple.trials.single.partner", "multiple.trials.multiples.partners", "one.try.each.partner" and "multiple.trials.each.partner". See section Methods for details  (default method.rewiring = "one.try.single.partner").
one.second.extinct.mod <- function(web, participant = "higher", method = "abun", ext.row = NULL, ext.col = NULL, 
                                   rewiring = FALSE, probabilities.rewiring1 = NULL, probabilities.rewiring2 = NULL,
                                   method.rewiring = "one.try.single.partner") {
  dead <- matrix(nrow = 0, ncol = 3)
  colnames(dead) <- c("no", "ext.lower", "ext.higher")
  m2 <- web
  i <- 1
  METHOD.REWIRING = c("one.try.single.partner", "multiple.trials.single.partner", "multiple.trials.multiples.partners", "one.try.each.partner", "multiple.trials.each.partner")
  method.rewiring <- pmatch(method.rewiring, METHOD.REWIRING)
  if (length(method.rewiring) > 1) {
    stop("\n Only one argument is accepted in method.rewiring \n")
  }
  if (is.na(method.rewiring)) {
    stop("\n Invalid method.rewiring \n")
  }
  if(method.rewiring == 4 | method.rewiring == 5){
    keep.trying <- TRUE
  } else {
    keep.trying <- FALSE
  }
  method.rewiring <- ifelse(method.rewiring == 4, 1, ifelse(method.rewiring == 5, 2, method.rewiring))
  if(rewiring){
    if(any(web%%1!=0)){
      stop("\n If rewiring is TRUE the web must must contain only integers \n")
    }
    if(is.null(rownames(web)) | is.null(colnames(web))){
      stop("\n If rewiring is TRUE the web must must rownames and colnames\n")
    }
    if(is.null(probabilities.rewiring1) | is.null(probabilities.rewiring2)){
      stop("\n If rewiring is TRUE probabilities.rewiring1 and probabilities.rewiring1 must not be NULL\n")
    }
  }
  repeat {
    ext.temp <- extinction.mod(m2, participant = participant, method = method, ext.row = ext.row, ext.col = ext.col)
    if(rewiring){
      if(!is.null(ext.temp$rexcl)){
        sp.ext <- rownames(ext.temp$rexcl)
        sp.try.rewiring <- which(ext.temp$rexcl>0)
        sp.surv <- seq_len(nrow(ext.temp$web))
        sp.surv <- sp.surv[-1*which(rownames(ext.temp$web) %in% sp.ext)]
        for(jj in sp.try.rewiring){
          sp.surv.temp <- sp.surv
          go <- TRUE
          m <- 0
          if(method.rewiring == 1 | method.rewiring == 3){
            trials <- 1
          } else {
            trials <- ext.temp$rexcl[1, jj]
          }
          while (go) {
            m <- m+1
            sp.surv.prob1 <- probabilities.rewiring1[sp.surv.temp, jj]
            sp.add <- sample(as.character(sp.surv.temp), 1, prob = sp.surv.prob1)
            sp.surv.prob2 <- probabilities.rewiring2[as.numeric(sp.add), jj]
            n.add <- rbinom(1, trials, sp.surv.prob2)
            if(n.add>0){
              ext.temp$web[as.numeric(sp.add), jj] <- ext.temp$web[as.numeric(sp.add), jj]+n.add
            }
            if(method.rewiring == 1 | method.rewiring == 2){
              if(!keep.trying){
                go <- FALSE  
              } else{
                if((method.rewiring == 1 & n.add>0) | (method.rewiring == 2 & n.add == trials)){
                  go <- FALSE
                } else{
                  sp.surv.temp <- sp.surv.temp[-1*which(sp.surv.temp%in% sp.add)]
                  trials <- trials-n.add
                  if(length(sp.surv.temp)<1){
                    go <- FALSE
                  }
                }
              }
            } else {
              if(ext.temp$rexcl[1, jj] == m){
                go <- FALSE
              }
            }
          }
        }
      }
      if(!is.null(ext.temp$cexcl)){
        sp.ext <- colnames(ext.temp$cexcl)
        sp.try.rewiring <- which(ext.temp$cexcl>0)
        sp.surv <- seq_len(ncol(ext.temp$web))
        sp.surv <- sp.surv[-1*which(colnames(ext.temp$web) %in% sp.ext)]
        for(ii in sp.try.rewiring){
          sp.surv.temp <- sp.surv
          go <- TRUE
          m <- 0
          if(method.rewiring == 1 | method.rewiring == 3){
            trials <- 1
          } else {
            trials <- ext.temp$cexcl[ii, 1]
          }
          while (go) {
            m <- m+1
            sp.surv.prob1 <- probabilities.rewiring1[ii, sp.surv.temp]
            sp.add <- sample(as.character(sp.surv.temp), 1, prob = sp.surv.prob1)
            sp.surv.prob2 <- probabilities.rewiring2[ii, as.numeric(sp.add)]
            n.add <- rbinom(1, trials, sp.surv.prob2)
            if(n.add>0){
              ext.temp$web[ii, as.numeric(sp.add)] <- ext.temp$web[ii, as.numeric(sp.add)]+n.add
            }
            if(method.rewiring == 1 | method.rewiring == 2){
              if(!keep.trying){
                go <- FALSE  
              } else{
                if((method.rewiring == 1 & n.add>0) | (method.rewiring == 2 & n.add == trials)){
                  go <- FALSE
                } else{
                  sp.surv.temp <- sp.surv.temp[-1*which(sp.surv.temp%in% sp.add)]
                  trials <- trials-n.add
                  if(length(sp.surv.temp)<1){
                    go <- FALSE
                  }
                }
              }
            } else {
              if(ext.temp$cexcl[ii, 1] == m){
                go <- FALSE
              }
            }
          }
        }
      }
    }
    n <- ext.temp$web
    dead <- rbind(dead, c(i, attributes(m2 <- empty(n, count = TRUE))$empty))
    if (participant == "lower" & NROW(m2) < 2) 
      break
    if (participant == "higher" & NCOL(m2) < 2) 
      break
    if (participant == "both" & min(dim(m2)) < 2) 
      break
    if (any(dim(n) == 1)) 
      break
    if (method == "external") {
      ext.col[ext.col > ext.col[1]] <- ext.col[ext.col > ext.col[1]] - 1
      ext.row[ext.row > ext.row[1]] <- ext.row[ext.row > ext.row[1]] - 1
      ext.row <- ext.row[-1]
      ext.col <- ext.col[-1]
    }
    i <- i + 1
  }
  dead2 <- rbind(dead, c(NROW(dead) + 1, NROW(m2), NCOL(m2)))
  if (participant == "lower" & method == "degree") {
    if (length(table(dead[, 2])) > 1) 
      dead2[, 2] <- 1
  }
  if (nrow(dead) + 1 != nrow(dead2)) 
    stop("PANIC! Something went wrong with the extinct sequence! Please contact the author to fix this!!")
  if (participant == "lower") 
    supposed.length <- NROW(web)
  if (participant == "higher") 
    supposed.length <- NCOL(web)
  if (participant == "both") 
    supposed.length <- NROW(dead2)
  if (NROW(dead2) != supposed.length) {
    missing <- supposed.length - NROW(dead2)
    addit1 <- (NROW(dead2) + 1):(NROW(dead2) + missing)
    addit2n3 <- rep(0, times = missing)
    dead2 <- rbind(dead2, as.matrix(data.frame(addit1, addit2n3, addit2n3)))
  }
  out <- dead2
  class(out) <- "bipartite"
  attr(out, "exterminated") <- c("both", "lower", "higher")[pmatch(participant, c("both", "lower", "higher"))]
  attr(out, "exterminated")
  return(out)
}

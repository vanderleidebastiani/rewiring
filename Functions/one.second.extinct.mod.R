# Small change in internal function "one.second.extinct" of bipartite package
# The difference is that allow rewiring the interation after each step of extinction
# New arguments:
# rewiring - Allow rewiring (default rewiring = FALSE)
# probabilities.rewiring1 - A matrix with probabilities of rewiring, must be the same dimensions of web (default probabilities.rewiring1 = NULL). See step ii?
# probabilities.rewiring2 - A matrix with probabilities of rewiring, must be the same dimensions of web (default probabilities.rewiring2 = NULL). See step iii?
# method.rewiring = Type of method used to trial rewiring, partial match to "one.try.single.interaction", "one.try.n.interactions", "multiple.trials" (Default method.rewiring = "one.try.n.interactions"). See text?
one.second.extinct.mod <- function(web, participant = "higher", method = "abun", ext.row = NULL, ext.col = NULL, 
                                   rewiring = FALSE, probabilities.rewiring1 = NULL, probabilities.rewiring2 = NULL,
                                   method.rewiring = "one.try.n.interactions") {
  dead <- matrix(nrow = 0, ncol = 3)
  colnames(dead) <- c("no", "ext.lower", "ext.higher")
  m2 <- web
  i <- 1
  METHOD.REWIRING = c("one.try.single.interaction", "one.try.n.interactions", "multiple.trials")
  method.rewiring <- pmatch(method.rewiring, METHOD.REWIRING)
  if (length(method.rewiring) > 1) {
    stop("\n Only one argument is accepted in method.rewiring \n")
  }
  if (is.na(method.rewiring)) {
    stop("\n Invalid method.rewiring \n")
  }
  if(rewiring){
    if(any(web%%1!=0)){
      stop("\n If rewiring is TRUE the web must must contain only integers \n")
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
          go <- TRUE
          m <- 0
          sp.surv.prob1 <- probabilities.rewiring1[sp.surv, jj]
          if(method.rewiring == 1 | method.rewiring == 3){
            trials <- 1
          } else {
            trials <- ext.temp$rexcl[1, jj]
          }
          while (go) {
            m <- m+1
            sp.add <- sample(as.character(sp.surv), 1, prob = sp.surv.prob1)
            sp.surv.prob2 <- probabilities.rewiring2[as.numeric(sp.add), jj]
            n.add <- rbinom(1, trials, sp.surv.prob2)
            if(n.add>0){
              ext.temp$web[as.numeric(sp.add), jj] <- ext.temp$web[as.numeric(sp.add), jj]+n.add
            }
            if(method.rewiring == 1 | method.rewiring == 2){
              go <- FALSE
            } else {
              if(trials == m){
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
          go <- TRUE
          m <- 0
          sp.surv.prob1 <- probabilities.rewiring1[ii, sp.surv]
          if(method.rewiring == 1 | method.rewiring == 3){
            trials <- 1
          } else {
            trials <- ext.temp$cexcl[ii, 1]
          }
          while (go) {
            m <- m+1
            sp.add <- sample(as.character(sp.surv), 1, prob = sp.surv.prob1)
            sp.surv.prob2 <- probabilities.rewiring2[ii, as.numeric(sp.add)]
            n.add <- rbinom(1, trials, sp.surv.prob2)
            if(n.add>0){
              ext.temp$web[ii, as.numeric(sp.add)] <- ext.temp$web[ii, as.numeric(sp.add)]+n.add
            }
            if(method.rewiring == 1 | method.rewiring == 2){
              go <- FALSE
            } else {
              if(trials == m){
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
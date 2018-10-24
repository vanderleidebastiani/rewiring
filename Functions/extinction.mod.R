# Small change in internal function "extinction" of bipartite package
# The diference is that return rexcl and cexcl (default NULL to both) with complele row (rexcl) or column (cexcl) extinct step by step
extinction.mod <- function (web, participant = "both", method = "random", ext.row = NULL, ext.col = NULL) 
{
  RES <- list(rexcl = NULL, cexcl = NULL)
	partis <- c("lower", "higher", "both")
	partis.match <- pmatch(participant, partis)
	if (is.na(partis.match)) 
		stop("Choose participant: lower/higher/both.\n")
	meths <- c("random", "abundance", "degree", "external")
	meths.match <- pmatch(method, meths)
	if (is.na(meths.match)) 
		stop("Choose extinction method: random/abundance/degree.\n")
	nr <- NROW(web)
	nc <- NCOL(web)
	if (partis.match == 3 & meths.match == 1) 
		partis.match <- sample(2, 1)
	if (meths.match == 1) {
		rexcl <- sample(nr, 1)
		cexcl <- sample(nc, 1)
		if (partis.match == 1){
		  RES$rexcl <- web[rexcl, , drop = FALSE]
			web[rexcl, ] <- 0
		}
		if (partis.match == 2) {
		  RES$cexcl <- web[, cexcl, drop = FALSE]
			web[, cexcl] <- 0
		}
	}
	if (meths.match == 2) {
		web <- web[sample(1:nrow(web)), sample(1:ncol(web)), drop = FALSE]
		rseq <- order(rowSums(web))
		cseq <- order(colSums(web))
		if (partis.match == 1){
		  RES$rexcl <- web[rseq[1], , drop = FALSE]
			web[rseq[1], ] <- 0
		}
		if (partis.match == 2){ 
		  RES$cexcl <- web[, cseq[1], drop = FALSE]
			web[, cseq[1]] <- 0
		}
		if (partis.match == 3) {
			if (min(rowSums(web)) < min(colSums(web))) {
			  RES$rexcl <- web[rseq[1], , drop = FALSE]
				web[rseq[1], ] <- 0
			}
			else {
				if (min(rowSums(web)) > min(colSums(web))) {
				  RES$cexcl <- web[, cseq[1], drop = FALSE]
					web[, cseq[1]] <- 0
				}
				else {
					if (sample(2, 1) == 1) {
					  RES$rexcl <- web[rseq[1], , drop = FALSE]
						web[rseq[1], ] <- 0
					}
					else {
					  RES$cexcl <- web[, cseq[1], drop = FALSE]
						web[, cseq[1]] <- 0
					}
				}
			}
		}
	}
	if (meths.match == 3) {
		if (partis.match == 1) {
			sequ <- rowSums(web > 0)
			which.ex <- which(sequ == max(sequ))
			if (length(which.ex) > 1) {
				ex <- sample(which.ex, size = 1)
			}
			else {
				ex <- which.ex
			}
			RES$rexcl <- web[ex, , drop = FALSE]
			web[ex, ] <- 0
		}
		if (partis.match == 2) {
			sequ <- colSums(web > 0)
			which.ex <- which(sequ == max(sequ))
			if (length(which.ex) > 1) 
				ex <- sample(which.ex, size = 1)
			else ex <- which.ex
			RES$cexcl <- web[, ex, drop = FALSE]
			web[, ex] <- 0
		}
	}
	if (meths.match == 4) {
		rseq <- ext.row
		cseq <- ext.col
		if (partis.match == 1){ 
		  RES$rexcl <- web[rseq[1], , drop = FALSE]
			web[rseq[1], ] <- 0
		}
		if (partis.match == 2) {
		  RES$cexcl <- web[, cseq[1], drop = FALSE]
			web[, cseq[1]] <- 0
		}
	}
	RES$web <- web
	return(RES)
}
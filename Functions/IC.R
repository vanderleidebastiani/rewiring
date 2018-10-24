# Small function to facilitate the calc of 95 percent confidence interval
IC <- function(X){
	a <- mean(X)
	s <- sd(X)
	n <- length(X)
	error <- qt(0.975, df = n-1)*s/sqrt(n)
	res <- c(mean = a, lower = a-error, upper = a+error)
	return(res)
}
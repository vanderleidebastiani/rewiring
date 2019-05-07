# Function to facilitate the calculation confidence interval (default conf.level = 0.95)
IC <- function(X, conf.level = 0.95){
  m <- mean(X)
  s <- sd(X)
  n <- length(X)
  se <- s/sqrt(n)
  error <- qt((1-conf.level)/2, df = n-1, lower.tail = FALSE)*se
  res <- c(mean = m, sd = s, lower = m-error, upper = m+error)
  return(res)
}

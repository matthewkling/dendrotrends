

logit <- function(x) log(x / (1-x))

inv_logit <- function(x) 1 / (1 + exp(-x))

decile <- function(x) floor(rank(x) / (length(x) + 1) * 10)

project <- function(x, y, p){
      mx <- mean(x)
      x <- x - mx
      beta <- sum(1/sum(x^2) * x * y) # formula for regression coef if mean(x) == 0
      return(mean(y) + beta * (p - mx))

      # demonstrate that lm code is legit
      # x <- rnorm(3)
      # x <- x - mean(x)
      # y <- rnorm(3, 10)
      # coef(lm(y ~ x))
      # sum(1/sum(x^2) * x * y)
      # x <- matrix(x, ncol = 1)
      # solve(t(x) %*% x) %*% t(x) %*% y
}



prj_logit <- function(x, y, x2, offset= 1/10000){
      if(length(unique(x)) == 1) return(rep(mean(y), length(x2)))
      y <- y * (1 - offset) + offset/2
      inv_logit(project(x, logit(y), x2))
}


ann2multi <- function(x, t) 1 - (1 - x) ^ t

multi2ann <- function(x, t) 1 - ((1 - x) ^ (1/t))

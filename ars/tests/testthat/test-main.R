main <- function(N){
  check_pval <- function(output, fun.name){
    if(output$p.value <= 0.05){
      print(paste("did NOT create a good sample with", fun.name,
                  ", pvalue=", round(output$p.value, 3)))
    } else {
      print(paste("passed with", fun.name))
    }
  }
  
  test_dist <- function(FUN, fun.name, a, b, N, rFUN, ...){
    out <- try(ars(N=N, f=FUN, a=a, b=b),
               silent=TRUE)
    if(inherits(out, "try-error")) {
      print(paste("fails to evaluate with", fun.name))
    } else {
      y1 <- rFUN(N,...)
      test <- ks.test(out$sample, y1)
      check_pval(test, fun.name)
    }
  }
  
  ## *******************************************************
  ## i. Standard Normal Distribution
  fun.name <- "normal"
  dnor <- function(x){
    return((1/(sqrt(2*pi)))*exp(-(x^2)/2))
  }
  test_dist(FUN=dnor,
            fun.name=fun.name,
            a=-10, b=10, N=N, rFUN=rnorm)
  ## *******************************************************
  ## ii. Laplace Distribution
  fun.name <- "laplace"
  dlaplace <- function(x, m = 0, s = 1){
    return(exp(-abs(x-m)/s)/(2*s))
  }
  rlaplace <- function(N){
    U <- runif(N, -0.5, 0.5)
    y <- sign(U)*rexp(N, 1)
    return(y)
  }
  
  test_dist(FUN=dlaplace,
            fun.name=fun.name,
            a=-10, b=10, N=N, rFUN=rlaplace)
  ## *******************************************************
  ## iii. Logistic Distribution
  fun.name <- "logistic"
  dlogistic <- function(x){
    return(exp(x)/(1+exp(x))^2)
  }
  test_dist(FUN=dlogistic,
            fun.name=fun.name,
            a=-10, b=10, N=N, rFUN=rlogis)
  
  ## *******************************************************
  ## iv. Gamma Distribution (chi square)
  fun.name <- "gamma"
  dgam <- function(x, theta=2, k=2){
    return((1/(gamma(k)*theta^k))*(x^(k-1))*exp(-x/theta))
  }
  test_dist(FUN=dgam,
            fun.name=fun.name,
            a=0.001, b=1000, N=N, rFUN=rgamma,
            shape = 2, scale = 2)
  
  ## *******************************************************
  ## v.  uniform Distribution (chi square)
  fun.name <- "uniform"
  dun <- function(x){
    return(dunif(x, min=0, max=1))
  }
  test_dist(FUN=dun,
            fun.name=fun.name,
            a=0, b=1, N=N, rFUN=runif)
  
  ## *******************************************************
  ## vi.  beta
  fun.name <- "beta"
  dbet <- function(x){
    return(dbeta(x, shape1 = 3, shape2 = 2, ncp = 0, log = FALSE))
  }
  test_dist(FUN=dbet,
            fun.name=fun.name,
            a=0.01, b=0.99, N=N, rFUN=rbeta,
            shape1=3, shape2=2)
  
  ## *******************************************************
  ## vii. non-concave distribution
  fun.name <- "non-log concave exp(x^2)"
  de <- function(x){
    return(exp(x^2))
  }
  test_dist(FUN=de,
            fun.name=fun.name,
            a=-1, b=1, N=N, rFUN=NA)
}

main(1000)
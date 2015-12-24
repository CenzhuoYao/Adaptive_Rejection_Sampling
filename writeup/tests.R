setwd('~/GitHub/243project/writeup')
source('~/GitHub/243project/ars/R/tmp_debug.R')

## nice little pdf function
pdf_f <- function(f, file, ...) {
  cat(sprintf("Writing %s\n", file))
  pdf(file, ...)
  on.exit(dev.off())
  f()
}

plot_densities <- function(){
  layout(matrix(1:6, nrow=3))
  par(oma=c(2, 1, 1, 1),
      mar=c(2, 2, 2, 1),
      mgp=c(2, 1, 0))
  ## *******************************************************
  ## i. Standard Normal Distribution
  x <- seq(-5, 5, len = N)
  dnor <- function(x){
    return((1/(sqrt(2*pi)))*exp(-(x^2)/2))
  }

  out.norm <- calc_sample(N=N, f=dnor, a=-5, b=5)
  plot(density(out.norm$sample), type="l", col="dodgerblue", lwd=2,
       main="Normal distribution", xlab="")
  y <- dnorm(x)
  points(x, y, type = "l", lwd=2)
  ## *******************************************************
  ## ii. Laplace Distribution
  x <- seq(-5, 5, len = N)
  dlaplace <- function(x, m = 0, s = 1){
    return(exp(-abs(x-m)/s)/(2*s))
  }

  out.lap <- calc_sample(N=N, f= dlaplace, a=-5, b=5)
  plot(density(out.lap$sample), type="l", col="dodgerblue", lwd=2,
       main="Laplace distribution", xlab="")
  y <- dlaplace(x, m=0, s=1)
  points(x, y, type = "l")
  ## *******************************************************
  ## iii. Logistic Distribution
  x <- seq(-10, 10, len = N)
  dlogistic <- function(x){
    return(exp(x)/(1+exp(x))^2)
  }

  out.log <- calc_sample(N=N, f=dlogistic, a=-10, b=10)
  plot(density(out.log$sample), type="l", col="dodgerblue", lwd=2,
       main="Logistic distribution", xlab="")
  y <- dlogistic(x)
  points(x, y, type = "l")

  ## *******************************************************
  ## iv. Gamma Distribution (chi square)
  x <- seq(0, 100, len = N)
  dgam <- function(x, theta=2, k=2){
    return((1/(gamma(k)*theta^k))*(x^(k-1))*exp(-x/theta))
  }

  out.gam <- calc_sample(N=N, f=dgam, a=0.001, b=1000)

  plot(density(out.gam$sample), type="l", col="dodgerblue", lwd=2,
       main="Gamma distribution", xlab="")
  y <- dgamma(x, shape = 2, scale = 2)
  points(x, y, type = "l")

  legend("topright", legend=c("Sample", "Theoretical"),
         col=c("dodgerblue", "black"), lwd=2, bty="n")

  ## *******************************************************
  ## v. Uniform Distribution
  dun <- function(x){
    return(dunif(x, min=0, max=1))
  }
  out.unif <- calc_sample(N=N, f=dun, a=0, b=1)
  plot(density(out.unif$sample), type="l", col="dodgerblue", lwd=2,
       main="Uniform distribution", xlab="")

  y <- dunif(x)
  points(x, y, type = "l")

  ## *******************************************************
  ## vi. Beta Distribution
  x <- seq(0, 1, len = N)
  
  dbet <- function(x, alpha=2, b=2){
    return((x^(alpha-1)*(1-x)^(b-1))/beta(alpha, b))
  }

  out.beta <- calc_sample(N=N, f=dbet, a=0.01, b=0.99)

  plot(density(out.beta$sample), type="l",
       col="dodgerblue", lwd=2,
       main="Beta distribution", xlab="")

  y <- dbet(x)
  points(x, y, type = "l")

}

N <- 10000
pdf_f(plot_densities, file="figures/densities.pdf",
      height=6, width=4)

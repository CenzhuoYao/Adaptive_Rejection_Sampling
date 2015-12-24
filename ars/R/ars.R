#' Adapative Rejection Sampling
#'
#' Adapative rejection sampling from log-concave functions. See Gilks and Wild (1992) for more details.
#' @param N sample size (or number of observations).
#' @param f function f(x) which is log-concave and is proporitional to the density you wish to sample from
#' @param a first starting point for which f(x) is defined. The default is Inf.
#' @param b second starting pointfor which f(x) is defined. The default is -Inf.
#' @param DFUN the derivative of f(x)
#' @param step value of bandwidth around est_mode used to initialize function. The default is 0.5.
#' @param est_mod estimated mode of f(x) after log if f(x) is unbounded. The default is 0.
#' @return \code{ars} returns a list containing the following components:
#' \item{f}{f(x)}
#' \item{fprime}{the derivative of f(x)}
#' \item{sample}{a vector of length N containing sampled values}
#' \item{hit_rate}{a value giving the proportion of accepted sample values over all evaluated values}
#' @author 
#' Code developed by Xinyue Zhou, Cenzhuo Yao, Katherine C. Ullman, and Lauren Catherine Ponisio, based on Gilks and Wild (1992).
#' @references 
#' Gilks, W.R., P. Wild. (1992) Adaptive Rejection Sampling for Gibbs Sampling, Applied Statistics 41:337–348. 
#' @export
#' @examples 
#' # Example 1: Sample from Standard Normal Distribution
#' library(ars)
#'  N <- 1000
#'  x <- seq(-5, 5, len = N)
#'  dnor <- function(x) {
#'    return((1/(sqrt(2*pi)))*exp(-(x^2)/2))
#'  }
#'  out.norm <- ars(N=N, f=dnor, a=-5, b=5)
#'  plot(density(out.norm$sample), type="l", col="dodgerblue", lwd=2,
#'       main="Normal distribution", xlab="")
#'  y <- dnorm(x)
#'  points(x, y, type = "l", lwd=2)
#' # Example 2 : Sample from Laplace Distribution
#'  x <- seq(-5, 5, len = N)
#'  dlaplace <- function(x, m = 0, s = 1) {
#'  return(exp(-abs(x-m)/s)/(2*s))
#'  }
#' out.lap <- ars(N=N, f= dlaplace, a=-5, b=5)
#' plot(density(out.lap$sample), type="l", col="dodgerblue", lwd=2,
#'     main="Laplace distribution", xlab="")
#' y <- dlaplace(x, m=0, s=1)
#' points(x, y, type = "l")

## main function to implement adaptive rejection algorithm
## Input
## N - number: the number of desired samples
## f - function: the density function users interested in
## a, b - number: the optional starting values
## DFUN - function: the optional fprime
## step - number: value of bandwidth around est_mode used to initialize function, optional.
## est_mod - number: estimated mode of f(x) after log if f(x) is unbounded, optional.

## Output - list
## f - function: the density function users interested in
## fprime - function: fprime we finally used
## sample - vector: the samples drawn from the interested density function
## hit_rate - number: a value giving the proportion of accepted sample values over all evaluated values

ars <- function(N, f,
                a = -Inf,
                b = Inf,
                DFUN = NA,
                step = 0.5,
                est_mod = 0){
  ## check if it is unbounded or not 
  if(class(f) != "function"){
    stop('please provide f as an function', call. = FALSE)
  }
  if(!is.na(DFUN)){
    if(class(DFUN) != "function"){
      warning('DFUN is not an function, re-creating from FUN',
              call. = FALSE, immediate.= TRUE)
      DFUN <- NA
    }
  }
  if(a == b){
    stop('please provide different a and b', call. = FALSE)
  }
  
  ##log of density function: right now, users can input density function
  FUN <- function(x,fun = f){
    return(log(fun(x)))
  }
  ##differentiate it using the function we define
  Derv_final <- function(x,fun=FUN){
    return(Derv(x,fun,a,b))
  }
  
  if(!is.na(DFUN)){
    DD<- DFUN
  } else{
    DD <- Derv_final
  }
  
  
  ## initial return sample
  ret <- rep(NA, N)
  ##count step
  count <- 0 
  info <- matrix(NA, nrow = as.integer((N^(1/3) + 2)),
                 ncol = 3)
  z <- c(NA, rep(NA, as.integer((N^(1/3) + 2))))
  ## matrix for x, f, fprime values, initialize at the expected length
  ## from Gilks et al. 1992 with a little extra
  ## how to initialize the info matrix
  if (a != -Inf && b != Inf){
    info <- help_init_info(a, b, info, FUN, DD)
    z <- init_z(x_star, a, b, info, z)
    ## for uniform case:
    if ((abs(info[2, 3]) < 1e-6)&&(abs(info[1, 3])< 1e-6)){
      return (list(f = FUN,
                   fprime = DD,
                   sample = runif(n=N ,a, b),
                   info = info[1:2,]))
    }
  } 
  
  if (a == -Inf && b != Inf){
    ## if est_mod out of the support, then
    if (est_mod > b) {
      est_mod <- b - step
    }
    aa <- est_mod
    test <- DD(est_mod)
    while (-Inf < test && test <= 0 && count <=50){
      aa <- aa - step
      test <- DD(aa)
      count = count + 1
    }
    info <- help_init_info(aa, b, info, FUN, DD)
    z <- init_z(x_star, a, b, info, z)
  }
  
  
  if (a != -Inf && b == Inf){
    if (est_mod < a) { 
      est_mod = a + step
    }
    bb <- est_mod
    test <- DD(bb)
    while (0 <=  test && test < Inf && count <=50){
      bb <- bb + step
      test <- DD(bb)
      count = count + 1
    }
    info <- help_init_info(a, bb, info, FUN, DD)
    z <- init_z(x_star, a, b, info, z)
  }
  
  if (a == -Inf && b == Inf){
    aa <- est_mod - step 
    bb <- est_mod + step
    test1 <- DD(aa)
    test2 <- DD(bb)
    while (-Inf < test1 && test1 <= 0 && count <= 50 ){
      aa <- aa - step
      test1 <- DD(aa)
      count = count + 1
    }
    while (0 <=  test2 && test2 < Inf && count <= 100){
      bb <- bb + step
      test <- DD(bb)
      count = count + 1
    }
    info <- help_init_info(aa, bb, info, FUN, DD)
    z <- init_z(x_star, a, b, info, z)
  }
  
  if (count >= 50) {
    stop ("est_mod is not valid, initial point cannot be found. Try another est_mod",.call = FALSE)
  }
  
  ## starting index
  itt1 <- 3
  itt2 <- 4
  sample <- it <- 0
  num_sample <- 1
  
  ## sample step
  while(sample < N){
    it <- it + num_sample
    hit_rate <- 1 - (itt1 - 3)/it
    clean_info <- info[1:(itt1-1),]
    ## check log-concavity of f(x)
    if(!check_concave(clean_info)) {
      stop('Input is not log-concave function', call. = FALSE)
    }
    clean_z <- z[1:(itt2-1)]
    ##sample!!!
    if (hit_rate == 1) {
      num_sample <- 1
    } else {
      num_sample <- round(1/(1-hit_rate), digits = 0)}
    
    x <- sample_envelope(clean_info, clean_z, num_sample)
    ##print(x)
    ## check if x is defined on fx
    test_fx <- FUN(x)
    ## draw a new x if not
    if(sum(!is.finite(test_fx))) next
    ## check if this value of x have been drawn before. If that is the
    ## case, we can directly throw it into the info matrix or sample
    
    ## get the corresponding point on upper bound of lower bound
    ## whether we reject the sample
    w <- runif(num_sample)
    for (i in 1:num_sample){
      ind1 <- rev(which(clean_info[ ,1] < x[i]))[1]
      ind2 <- which(clean_info[ ,1] > x[i])[1]
      if ((is.na(ind1) || is.na(ind2))&& is.finite(a) && is.finite(b)) {
        next
      }
      
      ## calculate lower bound
      if (is.na(ind1) || is.na(ind2)) {
        lx <- -Inf
      } else {
        lx <- line_fun(clean_info[ind1, 1],
                       clean_info[ind1, 2],
                       clean_info[ind2, 1],
                       clean_info[ind2, 2], x[i])
      }
      if (is.na(ind1)){
        ux <- line_fun_p(clean_info[ind2, 1],
                         clean_info[ind2, 2],
                         clean_info[ind2, 3], x[i])
      } else {
        if(x[i] <= clean_z[ind1 + 1]){
          ux <- line_fun_p(clean_info[ind1, 1],
                           clean_info[ind1, 2],
                           clean_info[ind1, 3], x[i])
        } else{
          ux <- line_fun_p(clean_info[ind2, 1],
                           clean_info[ind2, 2],
                           clean_info[ind2, 3], x[i])
        } 
      }
      
      ## reject/accept sample
      if(w[i] <= exp(lx - ux)) {
        sample <- sample + 1
        ret[sample] <- x[i] 
        next
      } else{
        if(w[i] <= exp(test_fx[i] - ux)){
          sample <- sample + 1
          ret[sample] <- x[i]
        }
        info <- update_info(x[i], info, itt1, FUN, test_fx[i],
                            DD,
                            num_sample)
        itt1 <- itt1 + 1
        z <- update_z(x[i], z, itt2, info[1:(itt1-1), ])
        itt2 <- itt2 + 1
      } 
    }
  }
  ## return sample
  hit_rate <- length(ret)/it
  return(list(f=FUN,
              fprime=DD,
              sample = ret,
              hit_rate = 1 - (itt1 - 3)/it))
}
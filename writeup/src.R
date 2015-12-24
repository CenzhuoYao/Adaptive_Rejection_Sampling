rm(list=ls())
require(fields)

## @knitr Derv
## function for numerically approximating the derivative
## x is the point to evaluate the derivate at, FUN is the function,
## a and b are the bounds
## the function returns the approximation to the derivative at the point of
## evaluation
Derv <- function(x, FUN, a, b){
  if (x == a) {return ((FUN(x + 1e-8)-FUN(x))/1e-8)}
  if (x == b) {return ((FUN(x) - FUN(x - 1e-8))/1e-8)}
  if (a <= x && x <= b) {return((FUN(x + 1e-8)-FUN(x - 1e-8))/2e-8)}
}

## @knitr help_init_info
## function for initialize the info matrix in unbounded cases
## aa and bb are the bounds, info is the info matrix with x, fx, and
## f'x, FUN is the function of interest and DD is the derivative
## function
## the function returns the updated info matrix
help_init_info <- function (aa, bb, info, FUN, DD){
  info[1, 1] <- aa
  info[1, 2] <- FUN(info[1, 1])
  if(!is.finite(info[1, 2])) {
    stop('aa is not defined on f', call. = FALSE)
  }
  info[1, 3] <- DD(info[1, 1])
  info[2, 1] <- bb
  info[2, 2] <- FUN(info[2, 1])
  if(!is.finite(info[1, 2])) {
    stop('bb is not defined on f', call. = FALSE)
  }
  info [2, 3] <- DD(bb)
  return(info)
}

## @knitr update_info
## updates the info matrix with new x, f, fprime values, and orders
## them takes xstar is the proposed sample, info is the info matrix
## with x, fx, and f'x, itt1 is the index, FUN is the function of
## interact, test_fx is XXXXXX, DD is the derivative function,
## num_sample is XXXXXX
## returns the updated info matrix
update_info <- function(x_star, info, itt1, FUN,
                        test_fx,
                        DD,
                        num_sample){
  if(nrow(info) < itt1+num_sample+2) {
    tmp <- matrix(NA, nrow = nrow(info), ncol = 3)
    info <- rbind(info, tmp)
  }
  ind <- sum(info[1:(itt1-1), 1] < x_star)[1]
  info[(ind+2):itt1, ] = info[(ind+1):(itt1-1), ]
  info[(ind+1), ] = c(x_star,test_fx,DD(x_star))
  return(info)
}

## @knitr init_z
## calculate z values, and order them, info should be clean_info
init_z <- function (x_star, a, b, info, z ){
  calc_z_new <-function (info){
    return((info[2, 2] - info[1, 2] -info[2, 1]*info[2, 3] +
            info[1, 1]*info[1, 3])/(info[1, 3] - info[2, 3]))
  }
  z[1] = a ; z[3] = b
  if (abs(info[2,3]-info[1,3])<1e-6) {
    z[2] = (info[1,1] + info[2,1])/2
  } else {
    z[2] = calc_z_new(info)
  }
  return (z)
}

## @knitr update_z
## updates the z vector
update_z <- function(x_star, z, itt2, info){
  if (nrow(info) <= 2){
    stop("something wrong with info matrix", call. = FALSE)
  }
  
  if(length(z) < itt2) {
    z <- c(z, rep(NA, length(z)))
  }
  
  ind <- which(info[, 1] == x_star)[1]
  if((ind > 1) && (ind < nrow(info))){
    if (sum(abs(diff(info[(ind-1):(ind+1), 3]))<1e-6) ){
      z_new1 <- info[ind - 1,1] + (info[ind + 1,1]-info[ind - 1,1])/3
      z_new2 <- info[ind - 1,1] + (info[ind + 1,1]-info[ind - 1,1])*2/3
    } else {
      z_new1 <-(info[ind, 2] - info[ind - 1, 2] -
                info[ind, 1]*info[ind, 3] +
                info[ind - 1, 1]*info[ind - 1, 3])/
                  (info[ind - 1, 3] - info[ind, 3])
      z_new2 <-(info[ind, 2] - info[ind + 1, 2] -
                info[ind, 1]*info[ind, 3] +
                info[ind + 1, 1]*info[ind + 1, 3])/
                  (info[ind + 1, 3] - info[ind, 3])
    }
    
    if((round(z_new1,digits = 6) > round(z_new2,digits = 6))){
      stop('something wrong with updated z value, check log-convexity!',
           call. = FALSE)
    }
    z[(ind+2):itt2] = z[(ind+1):(itt2-1)]
    z[ind] = z_new1
    z[ind+1] = z_new2
  } 
  if (ind == 1){
    if (sum(abs(diff(info[(ind):(ind+1), 3]))<1e-6) ){
      z_new1 <- -Inf
      z_new2 <- info[ind,1] + (info[ind + 1,1]-info[ind,1])/2 #take a look?
    } else {
      z_new1 <- -Inf
      z_new2 <- (info[ind, 2] - info[ind + 1, 2] -
                 info[ind, 1]*info[ind, 3] +
                 info[ind + 1, 1]*info[ind + 1, 3])/
                   (info[ind + 1, 3] - info[ind, 3])
    }
    z[(ind+2):itt2] = z[(ind+1):(itt2-1)]
    z[ind] = z_new1
    z[ind+1] = z_new2
  }

  if (ind == nrow(info)){
    if (sum(abs(diff(info[(ind-1):ind, 3]))<1e-6) ){
      z_new1 <- info[ind - 1,1] + (info[ind,1]-info[ind - 1,1])/2
      z_new2 <- Inf #take a look?
    } else {
      z_new1 <-(info[ind, 2] - info[ind - 1, 2] -
                info[ind, 1]*info[ind, 3] +
                info[ind - 1, 1]*info[ind - 1, 3])/
                  (info[ind - 1, 3] - info[ind, 3])
      z_new2 <- Inf
    }
    z[(ind+2):itt2] = z[(ind+1):(itt2-1)]
    z[ind] = z_new1
    z[ind+1] = z_new2
  }
  return(z)
}

## @knitr check_concave
## checks log-convexity
## clean_info is the ordered info matrix
check_concave <- function(clean_info){
  if(is.null(clean_info[,3])) {return(TRUE)}
  return(prod(round(clean_info[,3][-1],5) <=
              round(clean_info[,3][-nrow(clean_info)], 5)))
}

## @knitr line_fun
## takes two points and returns the function of that line
line_fun <- function(x1, y1, x2, y2, x_star){
  if(x1 == x2) {stop('Two points in a vertical line')}
  y_star <- (y2 - y1)/(x2 - x1)*(x_star - x1) + y1
  return(y_star)
}

## @knitr line_fun_p
## takes give x and h'(x), returns the function of that line
line_fun_p <- function(x1, y1, hpx1, x_star){
  y_star <- hpx1*(x_star - x1) + y1
  return(y_star)
}

## @knitr sample_envelope
## returns a new sample given the info matrix and z
## input is clean info matrix and clean z
sample_envelope <- function(samp_info, samp_z, num_sample){
  ## segments of bounds
  p <- rep(NA, nrow(samp_info))
  p <- exp(samp_info[,2])*(exp((samp_z[-1] - samp_info[,1])*
                               samp_info[,3]) -
                           exp((samp_z[-nrow(samp_info) - 1] -
                                samp_info[,1])*
                               samp_info[,3]))/samp_info[,3]
  ## q for normalizing p
  q <- p/sum(p)
  ## which region we sample from
  w <- runif(num_sample)
  i <- vector(mode = 'numeric', length = num_sample)
  for (j in 1:num_sample){
    i[j] <- sum(cumsum(q) < w[j]) + 1
  }
  ## sample pi using inv cdf
  samp_x <- (log(p[i]*runif(num_sample)*samp_info[i, 3] +
                   exp(samp_info[i,2] +
                         (samp_z[i] - samp_info[i,1])*samp_info[i,3])) -
               samp_info[i,2])/samp_info[i,3] + samp_info[i,1]
  return(samp_x)
}

## @knitr main_function
## main function to implement adaptive rejection algorithm
## fp is the optional fprime, expressed an an expression
## N is the number of desired samples
## a and b are the starting values
calc_sample <- function(N, f,
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
      est_mod = b - step
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
    stop
    ("est_mod is not valid, initial point cannot be found. Try another est_mod",
     .call = FALSE)
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
      print(clean_info)
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
        info <- update_info(x[i], info, itt1, FUN, test_fx[i],DD, num_sample)
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
              hit_rate = hit_rate,
              info = info[1:(itt1-1),],
              z = z[1:(itt2-1)]))
}

## @knitr samp_test
sample_envelope <- function(samp_info, samp_z, num_sample){
  ## segments of bounds
  p <- rep(NA, nrow(samp_info))
  p <- exp(samp_info[,2])*(exp((samp_z[-1] - samp_info[,1])*
                               samp_info[,3]) -
                           exp((samp_z[-nrow(samp_info) - 1] -
                                samp_info[,1])*
                               samp_info[,3]))/samp_info[,3]
  ## q for normalizing p
  q <- p/sum(p)
  ## which region we sample from
  w <- runif(num_sample)
  i <- vector(mode = 'numeric', length = num_sample)
  for (j in 1:num_sample){
    i[j] <- sum(cumsum(q) < w[j]) + 1
  }
  ## sample pi using inv cdf
  samp_x <- (log(p[i]*runif(num_sample)*samp_info[i, 3] +
                   exp(samp_info[i,2] +
                         (samp_z[i] - samp_info[i,1])*samp_info[i,3])) -
               samp_info[i,2])/samp_info[i,3] + samp_info[i,1]
  return(samp_x)
}


## test using the exponential
set.seed(4)
a <- 0
b <- 1
lambda <- 1

samp_info <- matrix(c(a, b,
                      log(dexp(a, rate=lambda)),
                      log(dexp(b, rate= lambda)),
                      -lambda,
                      -lambda), nrow=2)
samp_z <- c(0, 0.5, Inf)

test.out <- sample_envelope(samp_info, samp_z, 10000)

out.test <- ks.test(test.out, pexp)

if(out.test$p.value > 0.05){
  print("sample matches exponential, sample_envelope unit test passed")
}else{
  print("sample does not match exponential, sample_envelope unit test failed")
}

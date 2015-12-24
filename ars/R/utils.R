
## function for numerically approximating the derivative
## x is the point to evaluate the derivate at, FUN is the function,
## a and b are the bounds
## the function returns the approximation to the derivative at the point of
## evaluation
Derv <- function(x, FUN, a, b){
  if (x==a) {return ((FUN(x + 1e-8)-FUN(x))/1e-8)}
  if (x==b) {return ((FUN(x)-FUN(x - 1e-8))/1e-8)}
  if (a <= x && x <= b) {return((FUN(x + 1e-8)-FUN(x - 1e-8))/2e-8)}
}

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

## updates the info matrix with new x, f, fprime values, and orders
## them takes xstar is the proposed sample, info is the info matrix
## with x, fx, and f'x, itt1 is the index, FUN is the function of
## interact, test_fx is fx, DD is the derivative function,
## num_sample is number of samples we draw at same time.
## returns the updated info matrix
update_info <- function(x_star, info, itt1,
                        FUN,
                        test_fx,
                        DD,
                        num_sample){
  if(nrow(info) < itt1 + num_sample + 2) {
    tmp <- matrix(NA, nrow = nrow(info), ncol = 3)
    info <- rbind(info, tmp)
  }
  ind <- sum(info[1:(itt1-1), 1] < x_star)[1]
  info[(ind+2):itt1, ] <- info[(ind+1):(itt1-1), ]
  info[(ind+1), ] <- c(x_star,test_fx,DD(x_star))
  return(info)
}

## calculate z values, and order them, info should be clean_info
init_z <- function(x_star, a, b, info, z ){
  calc_z_new <-function(info){
    return((info[2, 2] - info[1, 2] - info[2, 1]*info[2, 3] +
              info[1, 1]*info[1, 3])/(info[1, 3] - info[2, 3]))
  }
  z[1] <- a ; z[3] <- b
  if (abs(info[2,3]-info[1,3]) < 1e-6) {
    z[2] <- (info[1,1] + info[2,1])/2
  } else {
    z[2] <- calc_z_new(info)
  }
  return (z)
}

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
      z_new1 <- info[ind - 1, 1] + (info[ind + 1, 1] -
                                      info[ind - 1, 1])/3
      z_new2 <- info[ind - 1, 1] + (info[ind + 1, 1] -
                                      info[ind - 1, 1])*2/3
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
    
    if((round(z_new1, digits = 6) > round(z_new2, digits = 6))){
      stop('something wrong with updated z value, check log-convexity!',
           call. = FALSE)
    }
    z[(ind+2):itt2] <- z[(ind+1):(itt2-1)]
    z[ind] <- z_new1
    z[ind+1] <- z_new2
  } 
  if (ind == 1){
    if (sum(abs(diff(info[(ind):(ind+1), 3])) < 1e-6) ){
      z_new1 <- -Inf
      z_new2 <- info[ind,1] + (info[ind + 1, 1] - info[ind,1])/2 #take a look?
    } else {
      z_new1 <- -Inf
      z_new2 <- (info[ind, 2] - info[ind + 1, 2] -
                   info[ind, 1]*info[ind, 3] +
                   info[ind + 1, 1]*info[ind + 1, 3])/
        (info[ind + 1, 3] - info[ind, 3])
    }
    z[(ind+2):itt2] <- z[(ind+1):(itt2-1)]
    z[ind] <- z_new1
    z[ind+1] <- z_new2
  }
  
  if (ind == nrow(info)){
    if (sum(abs(diff(info[(ind-1):ind, 3])) < 1e-6) ){
      z_new1 <- info[ind - 1, 1] + (info[ind, 1] - info[ind - 1, 1])/2
      z_new2 <- Inf #take a look?
    } else {
      z_new1 <-(info[ind, 2] - info[ind - 1, 2] -
                  info[ind, 1]*info[ind, 3] +
                  info[ind - 1, 1]*info[ind - 1, 3])/
        (info[ind - 1, 3] - info[ind, 3])
      z_new2 <- Inf
    }
    z[(ind+2):itt2] <- z[(ind+1):(itt2-1)]
    z[ind] <- z_new1
    z[ind+1] - z_new2
  }
  return(z)
}

## checks log-convexity
## clean_info is the ordered info matrix
check_concave <- function(clean_info){
  if(is.null(clean_info[, 3])) {return(TRUE)}
  return(prod(round(clean_info[, 3][-1], 5) <=
                round(clean_info[, 3][-nrow(clean_info)], 5)))
}

## takes two points and returns the function of that line
line_fun <- function(x1, y1, x2, y2, x_star){
  if(x1 == x2) {stop('Two points in a vertical line')}
  y_star <- (y2 - y1)/(x2 - x1)*(x_star - x1) + y1
  return(y_star)
}

## takes give x and h'(x), returns the function of that line
line_fun_p <- function(x1, y1, hpx1, x_star){
  y_star <- hpx1*(x_star - x1) + y1
  return(y_star)
}

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

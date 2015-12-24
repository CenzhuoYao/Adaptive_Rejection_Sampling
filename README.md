# Adaptive_Rejection_Sampling
R Package: Adaptive Reject Sampling (`ars`)

The aim of `ars` is to quickly generate sample with hard-to-evaluate density function. The algorithms of adaptive rejection sampling is introduced with detail in Gilks et al. (1992). The method we implemented is tangent approach instead of secant approach. The package is called 'ars' and you can install and use in R.

## Download the latest version

You can track (and contribute to) development of `ars` at https://github.com/CenzhuoYao/Adaptive_Rejection_Sampling. To Install the release version of `ars`, download this package from github with `install("ars")`.

## Performance of sampling

We have tested several log-concave density function and found `ars` performed well for those distribution. The [write up](https://github.com/CenzhuoYao/Adaptive_Rejection_Sampling/blob/master/writeup/writeup.pdf) shows how we implement the algorithm and what improvements we made, as well as the performance of some regular density functions.

## Contributors

This package is written by four people: Lauren Ponisio, Katherine Ullman, Xinyue Zhou and Cenzhuo Yao. If you found any bug, please accept our apology and e-mail to Yao <cenzhuoyao@berkeley.edu> so we can fix it as soon as possible.

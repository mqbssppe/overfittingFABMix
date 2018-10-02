

R> library(RcppArmadillo)
R> RcppArmadillo.package.skeleton(name = "fabMix", code_files = "fabMix.R", example_code=FALSE)

T> cp fabMix.cpp fabMix/src/fabMix.cpp

meta edit description kai namespace

R> library(Rcpp)
R> compileAttributes("fabMix")

meta edit to man

meta R CMD buil check etc

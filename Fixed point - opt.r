
options(digits=20)
source("Functions.r")

# Results
ca <- 400
k <- 0.05
MAP <- 1825
pkx <- 0.5

f1 <- Vectorize(function(wLr){
  averBif1 <- Vectorize(function(wLi)averBif(wLi, wLr))
  optwLi <- optimize(averBif1, c(0.14, 0.206), tol=.Machine$double.eps, maximum=T)
  res <- abs(optwLi$maximum-wLr)
  return(res)
})

optwL <- optimize(f1, c(0.1778, 0.1782), tol=.Machine$double.eps)
optwL
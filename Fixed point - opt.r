
options(digits=20)
source("Functions.r")

# Results
ca <- 400
k <- 0.05
MAP <- 1825
pkx <- 0.5

f1 <- Vectorize(function(wLr){
  averBif1 <- Vectorize(function(wLi)averBif(wLi, wLr))
  res1 <- averBif1(wLr)
  res2 <- optimize(averBif1, c(0.1, 0.3), tol=.Machine$double.eps, maximum=T)$objective
  res <- abs(res2-res1)
  return(res)
})

optwL <- optimize(f1, c(0.176, 0.18), tol=.Machine$double.eps)#0.17773536484440103
optwL
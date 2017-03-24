
options(digits=20)
source("optwL/optwL - Functions.r")

# Results
ca <- 400
k <- 0.05
MAP <- 1825
pkx <- 0.25
c <- 2.64

h3 <- 25
optwL <- vector("list")

for(i in 1:10){
  d <- i
  optwL[[i]] <- try(uniroot(optwLif, c(0.1, 0.3), tol=.Machine$double.eps^0.25))
}
optwL

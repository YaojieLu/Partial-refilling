
options(digits=20)
source("Functions.r")

# Results
ca <- 400
c <- 2.64

k <- 0.05
MAP <- 1825
h3 <- 25
pkx <- 0.5

optwL <- vector("list")
for(i in 1:10){
  d <- i
  optwL[[i]] <- try(uniroot(optwLif, c(0.1, 0.3), tol=.Machine$double.eps^0.25))
}
optwL

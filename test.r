
options(digits=20)
source("Functions.r")

# Results
ca <- 400
c <- 2.64

k <- 0.1
MAP <- 912.5
h3 <- 1
pkx <- 0.25

optwL <- vector("list")
for(i in 15){
  d <- i
  optwL[[i]] <- try(uniroot(optwLif, c(0.1, 0.3), tol=.Machine$double.eps^0.25))
}
optwL

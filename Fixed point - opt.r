
options(digits=20)
source("Functions.r")

# Results
ca <- 400
k <- 0.05
MAP <- 1825
pkx <- 0.5

optwL <- uniroot(optwLif, c(0.177, 0.178), tol=.Machine$double.eps)
optwL
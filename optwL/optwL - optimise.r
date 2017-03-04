
options(digits=20)
source("optwL/optwL - Functions.r")

# Results
ca <- 400
k <- 0.05
MAP <- 1825
pkx <- 0.5

c <- 2.64
d <- 5

optwL <- uniroot(optwLif, c(0.156779197039201, 0.160349116248024), tol=.Machine$double.eps^0.25)
optwL

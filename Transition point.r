options(digits=20)
source("Functions.r")

# Results
ca <- 400
k <- 0.05
MAP <- 1825
pkx <- 0.5

wL <- 0.17
wLL <- wLLf(wL)
wLff <- ff(wL)

spf1 <- Vectorize(function(w)spf(w, wL))
uniroot(spf1, c(wLL, wLff), tol=.Machine$double.eps)

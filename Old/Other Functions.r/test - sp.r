
source("Functions.r")

# Initializing
ca <- 400
pkx <- 0.5

wL <- 0.1745
wLL <- wLLf(wL)
sp <- spf(wL)

f1 <- Vectorize(function(w)gswLf(w, wL))
f2 <- Vectorize(function(w)gsmaxfm(w, wL))

curve(f1, wL, 0.3)
abline(v=wLL, col="red")
abline(v=sp)
abline(h=f1(sp))
curve(f2, wL, 0.3, add=T, col="blue")

f3 <- Vectorize(spf)
curve(f3, 0.16, 0.19)
abline(a=0, b=1)


options(digits=20)
source("Functions.r")

# Results
ca <- 400
k <- 0.05
MAP <- 1825
pkx <- 0.5

#optwLf(0.15)
x <- seq(0.177, 0.179, by=0.0005)
res <- optwLf(x)

# Figure
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", mar=c(4, 4, 2, 2), mfrow=c(1, 1))
plot(x, res, xlim=c(min(x), max(x)), ylim=c(0.15, 0.2), type="p",
     xlab=expression(italic(w[Lr])), ylab=expression(Best~italic(w[Li])),
     cex.lab=1.3, col="blue", lwd=2)
abline(a=0, b=1)

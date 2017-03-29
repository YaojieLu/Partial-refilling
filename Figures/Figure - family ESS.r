
source("Functions.r")

# Initializing
ca <- 400
c <- 2.64
d <- 10
pkx <- 0.5
h3 <- 1

SA <- c(0.12681344)
wLL <- wLLf(SA)

# Figures
# ESS gs(w)
Cols <- rainbow(length(SA))
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(3.5, 4, 1, 1), mfrow=c(1, 1))
plot(0, 0, xlim=c(0, 1), ylim=c(0, 0.3), cex.lab=1.5,
     type="n", xaxt="n", yaxt="n", xlab=NA, ylab=NA)

axis(1, xlim=c(0, 1), pos=-1, lwd=2)
mtext(expression(italic(w)), side=1, line=2, cex=1.3)
axis(2, ylim=c(0, 15), pos=0, lwd=2)
mtext(expression(italic(g[s])~(mol~m^-2~s^-1)), side=2, line=1.9, cex=1.3)

for(i in 1:length(SA)){
  with(list(wL=SA[i], Col=Cols[i], wLL=wLL[i]), {
    f1 <- Vectorize(function(w)gswLf1(w, wL))
    curve(f1, wLL, 1, add=T, col=Col)
  })
}

# ESS px(ps)
Cols <- rainbow(length(SA))
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(3.5, 4, 1, 1), mfrow=c(1, 1))
plot(0, 0, xlim=c(-10, 0), ylim=c(-10, 0), cex.lab=1.5,
     type="n", xaxt="n", yaxt="n", xlab=NA, ylab=NA)

axis(1, xlim=c(-10, 0), pos=-10, lwd=2)
mtext(expression(italic(psi[s])~(MPa)), side=1, line=2, cex=1.3)
axis(2, ylim=c(-10, 0), pos=-10, lwd=2)
mtext(expression(italic(psi[x])~(MPa)), side=2, line=1.9, cex=1.3)

for(i in 1:length(SA)){
  with(list(wL=SA[i], Col=Cols[i], wLL=wLL[i]), {
    f1 <- Vectorize(function(w)pxf(w, gswLf1(w, wL), wL))
    w <- seq(wLL, 1, by=(1-wLL)/100)
    
    x <- psf(w)
    y <- f1(w)
    points(x, y, col=Col)
  })
}

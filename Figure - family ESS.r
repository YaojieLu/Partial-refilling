
source("Functions.r")

# Initializing
ca <- 400
pkx <- 0.5

SA <- c(0.22, 0.21, 0.20, 0.15, 0.1)
wLL <- wLLf(SA)

# Figure
Cols <- rainbow(length(SA))
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(3.5, 4, 1, 1), mfrow=c(1, 1))
plot(0, 0, xlim=c(0, 1), ylim=c(0, 0.3), cex.lab=1.5,
     type="n", xaxt="n", yaxt="n", xlab=NA, ylab=NA)

axis(1, xlim=c(0, 1), pos=-1, lwd=2)
mtext(expression(italic(w)),side=1, line=2, cex=1.3)
axis(2, ylim=c(0, 15), pos=0, lwd=2)
mtext(expression(italic(g[s])~(mol~m^-2~s^-1)), side=2, line=1.9, cex=1.3)

# Family ESS
for(i in 1:length(SA)){
  with(list(wL=SA[i], Col=Cols[i], wLL=wLL[i]), {
    f1 <- Vectorize(function(w)gswLf(w, wL))
    curve(f1, wLL, 1, add=T, col=Col)
    #sp <- spf(wL)
    #abline(h=f1(sp))
    #abline(v=sp)
  })
}

#for(i in 1:length(SA)){
#  with(list(wL=SA[i]), {
#    f2 <- Vectorize(function(w)gsmaxfm(w, wL))
#    curve(f2, wL, 1, add=T, lty=2)
#  })
#}
#

source("Functions.r")
wLdata <- read.csv("Data/optwL.csv")
c <- 2.64
pe <- -1.58*10^-3
h3SA <- c(25, 100, 25, 25)
dSA <- c(1, 5, 10)
Psi50fd <- P50f(dSA)
wLSA <- c(as.numeric(subset(wLdata, k==0.05 & MAP==1825 & d==dSA[1] & h3==25 & pkx==0.5, select=c("optwL"))),
          as.numeric(subset(wLdata, k==0.05 & MAP==1825 & d==dSA[2] & h3==25 & pkx==0.5, select=c("optwL"))),
          as.numeric(subset(wLdata, k==0.05 & MAP==1825 & d==dSA[3] & h3==25 & pkx==0.5, select=c("optwL"))))

# Figure
Cols <- c("purple", "orange", "forestgreen")
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(3.5, 3.5, 1, 3.5), mfrow=c(1, 1))

# gs
plot(0, 0,
     type="n", xaxt="n", yaxt="n", xlab=NA, ylab=NA,
     xlim=c(-12, 0), ylim=c(0, 0.4), cex.lab=1.3)

axis(1, xlim=c(-12, 0), pos=0, lwd=2)
mtext(expression(psi[s]~(MPa)),side=1,line=2.3, cex=1.3)
axis(2, ylim=c(0, 0.4), pos=-12, lwd=2)
mtext(expression(italic(g[s])~(mol~m^-2~s^-1)),side=2,line=1.8, cex=1.3)

for(i in 1:length(dSA)){
  pkx <- 0.5
  h3 <- 25
  d <- dSA[i]
  wL <- wLSA[i]
  wLL <- wLLf(wL)
  psLL <- psf(wLL)
  f1 <- function(ps)gswLpsf(ps, wL)
  curve(f1, psLL, pe, col=Cols[i], add=T)
}

# PLC
par(new=TRUE)

plot(0, 0,
     type="n", xaxt="n", yaxt="n", xlab=NA, ylab=NA,
     xlim=c(-12, 0), ylim=c(0, 100), cex.lab=1.3)

axis(4, ylim=c(0, 100), pos=0, lwd=2)
mtext("PLC (%)", side=4,line=2.1, cex=1.3)

for(i in 1:length(dSA)){
  pkx <- 0.5
  h3 <- 25
  d <- dSA[i]
  wL <- wLSA[i]
  wLL <- wLLf(wL)
  psLL <- psf(wLL)
  f2 <- Vectorize(function(ps)ESSPLCpsf(ps, wL)*100)
  curve(f2, psLL, pe, col=Cols[i], lty=2, add=T)
}

legend("topleft", title=expression(psi[x50]), legend=round(Psi50fd, 1), pch=c(19), col=Cols)
box()

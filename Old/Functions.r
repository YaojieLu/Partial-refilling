
# ps(w)
psf <- function(w, pe=-1.58*10^-3, b=4.38)pe*w^(-b)

# the original PLC(px)
PLCf <- function(px)1-exp(-(-px/d)^c)

# modified PLC
PLCfm1 <- function(px, wL){
  pxL <- psf(wL)
  res <- ifelse(px>pxL, PLCf(pxL)-(PLCf(pxL)-PLCf(px))*pkx, PLCf(px))
  return(res)
}

# P50
P50f <- Vectorize(function(d){
  f1 <- function(px)exp(-(-px/d)^c)-0.5
  res <- uniroot(f1, c(-10, 0), tol=.Machine$double.eps)$root
  return(res)
})

# modified gsmax
gsmaxfm <- function(w, wL,
                    a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                    h2=l*LAI/nZ*p/1000, kxmax=5){
  # modified PLC
  PLCfm <- function(x)PLCf(pxL)-(PLCf(pxL)-PLCf(x))*pkx
  # modified xylem conductance function
  kxfm <- function(x)kxmax*(1-PLCfm(x))
  
  f1 <- function(x)(ps-x)*h2*kxfm(x)/(h*VPD)
  
  pxL <- psf(wL)
  ps <- psf(w)
  res <- ifelse(pxL<ps, optimize(f1, c(pxL, ps), tol=.Machine$double.eps, maximum=T)$objective, 0)
  return(res)
}

# modified pxmin
pxminfm <- function(w, wL, a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                    h2=l*LAI/nZ*p/1000, kxmax=5){
  # modified PLC
  PLCfm <- function(x)PLCf(pxL)-(PLCf(pxL)-PLCf(x))*pkx
  # modified xylem conductance function
  kxfm <- function(x)kxmax*(1-PLCfm(x))
  
  f1 <- function(x)(ps-x)*h2*kxfm(x)/(h*VPD)
  
  pxL <- psf(wL)
  ps <- psf(w)
  res <- ifelse(pxL<ps, optimize(f1, c(pxL, ps), tol=.Machine$double.eps, maximum=T)$maximum, 0)
  return(res)
}

## modified pxmin
#pxminfm <- function(w, wL,
#                    LAI=1, nZ=0.5, p=43200, l=1.8e-5, VPD=0.02,
#                    h2=l*LAI/nZ*p/1000, kxmax=5){
#  # modified PLC
#  PLCfm <- function(x)PLCf(pxL)-(PLCf(pxL)-PLCf(x))*pkx
#  # modified xylem conductance function
#  kxfm <- function(x)kxmax*(1-PLCfm(x))
#  
#  f1 <- function(x)(ps-x)*h2*kxfm(x)
#  
#  pxL <- psf(wL)
#  ps <- psf(w)
#  wgsmaxpxL <- wgsmaxpxLf(wL)
#  res <- ifelse(pxL<ps, ifelse(w>wgsmaxpxL, optimize(f1, c(pxL, ps), tol=.Machine$double.eps, maximum=T)$maximum, pxL), ps)
#  return(res)
#  
#}

# xylem water potential function
pxf <- function(w, gs, wL,
                a=1.6, LAI=1, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                h2=l*LAI/nZ*p/1000, kxmax=5){
  # modified PLC
  PLCfm <- function(x)PLCf(pxL)-(PLCf(pxL)-PLCf(x))*pkx
  # modified xylem conductance function
  kxfm <- function(x)kxmax*(1-PLCfm(x))
  
  f1 <- function(x)(ps-x)*h2*kxfm(x)-h*VPD*gs
  
  pxL <- psf(wL)
  ps <- psf(w)
  pxmin <- pxminfm(w, wL)
  message(w, " ", gs, " ", wL)
  res <- ifelse(pxmin<ps, uniroot(f1, c(pxmin, ps), tol=.Machine$double.eps)$root, ps)
  return(res)
}

# where gs reaches its max at px=pxL
wgsmaxpxLf <- function(wL,
                       a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                       h2=l*LAI/nZ*p/1000, kxmax=5){
  f1 <- function(w){
    ps <- psf(w)
    px <- pxL
    res <- -((h2*kxmax*(exp((-(px/d))^c)*(-1+pkx)*(-1+PLCmax)*px+pkx*(px+c*ps*(-(px/d))^c-c*px*(-(px/d))^c)))/(exp((-(px/d))^c)*px))
  }
  
  pxL <- psf(wL)
  PLCmax <- PLCf(pxL)
  x <- try(uniroot(f1, c(wL, 1), tol=.Machine$double.eps)$root, silent=TRUE)
  res <- ifelse(is.numeric(x), x, 1)
  return(res)
}

# Af(gs)
Af <- function(gs, ca=400, Vcmax=50, cp=30, Km=703, Rd=1, LAI=1)LAI*1/2*(Vcmax+(Km+ca)*gs-Rd-((Vcmax)^2+2*Vcmax*(Km-ca+2*cp)*gs+((ca+Km)*gs+Rd)^2-2*Rd*Vcmax)^(1/2))

# modified mf(w, gs)
mfm <- function(w, gs, wL){
  # modified PLC
  PLCfm <- function(x)PLCf(pxL)-(PLCf(pxL)-PLCf(x))*pkx
  
  pxL <- psf(wL)
  px <- pxf(w, gs, wL)
  res <- h3*(PLCfm(px)-PLCfm(0))
  return(res)
}

# modified B(w, gs)
Bfm <- function(w, gs, wL)Af(gs)-mfm(w, gs, wL)

# switch point
spf <- function(wL,
                ca=400, Vcmax=50, cp=30, Km=703, Rd=1, LAI=1,
                a=1.6, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                h2=l*LAI/nZ*p/1000, kxmax=5){
  f1 <- function(w){
    dAdgsf <- function(gs)(1/2)*LAI*(ca+Km+((-ca^2)*gs-gs*Km^2-Km*Rd-2*cp*Vcmax-Km*Vcmax+ca*(-2*gs*Km-Rd+Vcmax))/sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))
    f2 <- function(px)h3*c*exp(-(-px/d)^c)*(-px/d)^(c-1)/d*pkx*((exp((-(px/d))^c)*h*px*VPD)/(h2*kxmax*(exp((-(px/d))^c)*(-1+pkx)*(-1+PLCmax)*px+pkx*(px+c*ps*(-(px/d))^c-c*px*(-(px/d))^c))))
    
    ps <- psf(w)
    gsmax <- gsmaxfm(w, wL)
    pxmin <- pxminfm(w, wL)
    res <- dAdgsf(gsmax)-f2(pxmin)
    return(res)
  }
  
  pxL <- psf(wL)
  PLCmax <- PLCf(pxL)
  wgsmaxpxL <- wgsmaxpxLf(wL)
  x <- try(uniroot(f1, c(wL, wgsmaxpxL*(1-1e-10)), tol=.Machine$double.eps)$root, silent=TRUE)
  res <- ifelse(wgsmaxpxL>wL, ifelse(is.numeric(x), x, ifelse(f1(wL)>f1(wgsmaxpxL), 1, wL)), wL)
  return(res)
}

# family ESS 1
gswLf1 <- function(w, wL){
  Bfm1 <- function(gs)Bfm(w, gs, wL)
  gsmaxfm1 <- function(w)gsmaxfm(w, wL)
  
  res <- ifelse(0<gsmaxfm1(w), optimize(Bfm1, c(0, gsmaxfm1(w)), tol=.Machine$double.eps, maximum=T)$maximum, 0)
  return(res)
}

# family ESS
gswLf <- function(w, wL){
  Bfm1 <- function(gs)Bfm(w, gs, wL)
  gsmaxfm1 <- function(w)gsmaxfm(w, wL)
  
  sp <- spf(wL)
  res <- ifelse(w>sp, ifelse(0<gsmaxfm1(w), optimize(Bfm1, c(0, gsmaxfm1(w)), tol=.Machine$double.eps, maximum=T)$maximum, 0), gsmaxfm1(w))
  return(res)
}

# family ESS A(w)
AwLf <- function(w, wL)Af(gswLf1(w, wL))

# family ESS B(w)
BwLf <- function(w, wL)Bfm(w, gswLf(w, wL), wL)

# ESS g1(ps)
ESSg1psf <- Vectorize(function(ps, wL, VPD=0.02, a=1.6){
  f1 <- function(w)psf(w)-ps
  w <- uniroot(f1, c(0.001, 1), tol=.Machine$double.eps)$root
  res <- sqrt(VPD*100)*(ca*gswLf1(w, wL)/(a*AwLf(w, wL))-1)
  return(res)
})

# LHS endpoint
wLLf <- Vectorize(function(wL){
  BwLf1 <- Vectorize(function(w)BwLf(w, wL))
  res <- uniroot(BwLf1, c(wL, 1), tol=.Machine$double.eps)$root
  return(res)
})

# averB for invader
averBif <- function(wLi, wLr,
                    a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                    pe=-1.58*10^-3, b=4.38, h2=l*LAI/nZ*p/1000, kxmax=5,
                    gamma=1/((MAP/365/k)/1000)*nZ){
  wLLr <- wLLf(wLr)
  wLLi <- wLLf(wLi)
  spr <- spf(wLr)
  spi <- spf(wLi)
  
  gsmaxfmr <- function(w)gsmaxfm(w, wLr)
  gsmaxfmi <- function(w)gsmaxfm(w, wLi)
  
  gswLfr <- Vectorize(function(w)ifelse(w<wLLr, 0, ifelse(w>spr, gswLf1(w, wLr), gsmaxfmr(w))))
  gswLfi <- Vectorize(function(w)ifelse(w<wLLi, 0, ifelse(w>spi, gswLf1(w, wLi), gsmaxfmi(w))))
  
  Evf <- function(w)h*VPD*gswLfr(w)
  Lf <- function(w)Evf(w)+w/200
  rLf <- function(w)1/Lf(w)
  integralrLf <- Vectorize(function(w)integrate(rLf, w, 1, rel.tol=.Machine$double.eps^0.25)$value)
  fnoc <- function(w)1/Lf(w)*exp(-gamma*w-k*integralrLf(w))
  
  f1 <- Vectorize(function(w)Bfm(w, gswLfi(w), wLi)*fnoc(w))
  res <- integrate(f1, wLLi, 1, rel.tol=.Machine$double.eps^0.25)$value
  message(wLr, " ", wLi, " ", res)
  return(res)
}

# averages in monoculture
averBf <- function(wL,
                   a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                   pe=-1.58*10^-3, b=4.38, h2=l*LAI/nZ*p/1000, kxmax=5,
                   gamma=1/((MAP/365/k)/1000)*nZ){
  wLL <- wLLf(wL)
  sp <- spf(wL)
  
  gswLf <- Vectorize(function(w)ifelse(w<wLL, 0, ifelse(w>sp, gswLf1(w, wL), gsmaxfm(w, wL))))
  
  Evf <- function(w)h*VPD*gswLf(w)
  Lf <- function(w)Evf(w)+w/200
  rLf <- function(w)1/Lf(w)
  integralrLf <- Vectorize(function(w)integrate(rLf, w, 1, rel.tol=.Machine$double.eps^0.4)$value)
  fnoc <- function(w)1/Lf(w)*exp(-gamma*w-k*integralrLf(w))
  browser()
  res1 <- integrate(fnoc, 0, 1, rel.tol=.Machine$double.eps^0.4)#$value
  cPDF <- 1/res1$value
  
  fA <- Vectorize(function(w)Af(gswLf(w))*cPDF*fnoc(w))
  resA <- integrate(fA, wLL, 1, rel.tol=.Machine$double.eps^0.4)#$value
  fE <- Vectorize(function(w)Evf(w)*cPDF*fnoc(w))
  resE <- integrate(fE, wLL, 1, rel.tol=.Machine$double.eps^0.4)#$value
  res <- c(resA$value, resE$value*500*365)
  return(res)
}

optwLif <- Vectorize(function(wLr){
  averBif1 <- Vectorize(function(wLi)averBif(wLi, wLr))
  optwLi <- optimize(averBif1, c(0.11, 0.15), tol=.Machine$double.eps^0.25, maximum=T)
  res <- optwLi$maximum-wLr
  return(res)
})

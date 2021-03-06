
# ps(w)
psf <- function(w, pe=-1.58*10^-3, b=4.38)pe*w^(-b)

# the original PLC(px)
PLCf <- function(px, c=2.64, d=3.54)1-exp(-(-px/d)^c)

# modified gsmax
gsmaxfm <- function(w, wL,
                    a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                    h2=l*LAI/nZ*p/1000, kxmax=5, c=2.64, d=3.54){
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

# where gs reaches its max at px=pxL
wgsmaxpxLf <- function(wL,
                       a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                       h2=l*LAI/nZ*p/1000, kxmax=5, c=2.64, d=3.54){
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

# modified pxmin
pxminfm <- function(w, wL,
                    LAI=1, nZ=0.5, p=43200, l=1.8e-5, VPD=0.02,
                    h2=l*LAI/nZ*p/1000, kxmax=5){
  # modified PLC
  PLCfm <- function(x)PLCf(pxL)-(PLCf(pxL)-PLCf(x))*pkx
  # modified xylem conductance function
  kxfm <- function(x)kxmax*(1-PLCfm(x))
  
  f1 <- function(x)(ps-x)*h2*kxfm(x)
  
  pxL <- psf(wL)
  ps <- psf(w)
  wgsmaxpxL <- wgsmaxpxLf(wL)
  res <- ifelse(pxL<ps, ifelse(w>wgsmaxpxL, optimize(f1, c(pxL, ps), tol=.Machine$double.eps, maximum=T)$maximum, pxL), ps)
  return(res)
  
}

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
  res <- ifelse(pxmin<ps, uniroot(f1, c(pxmin, ps), tol=.Machine$double.eps)$root, ps)
  return(res)
}

# Af(gs)
Af <- function(gs, ca=400, Vcmax=50, cp=30, Km=703, Rd=1, LAI=1)LAI*1/2*(Vcmax+(Km+ca)*gs-Rd-((Vcmax)^2+2*Vcmax*(Km-ca+2*cp)*gs+((ca+Km)*gs+Rd)^2-2*Rd*Vcmax)^(1/2))

# modified mf(w, gs)
mfm <- function(w, gs, wL, h3=10){
  # modified PLC
  PLCfm <- function(x)PLCf(pxL)-(PLCf(pxL)-PLCf(x))*pkx
  
  pxL <- psf(wL)
  px <- pxf(w, gs, wL)
  res <- h3*(PLCfm(px)-PLCfm(0))
  return(res)
}

# modified B(w, gs)
Bfm <- function(w, gs, wL)Af(gs)-mfm(w, gs, wL)

# family ESS
gswLf <- function(w, wL){
  Bfm1 <- function(gs)Bfm(w, gs, wL)
  gsmaxfm1 <- function(w)gsmaxfm(w, wL)
  
  res <- ifelse(0<gsmaxfm1(w), optimize(Bfm1, c(0, gsmaxfm1(w)), tol=.Machine$double.eps, maximum=T)$maximum, 0)
  return(res)
}

# family ESS B(w)
BwLf <- function(w, wL)Bfm(w, gswLf(w, wL), wL)

# LHS endpoint
wLLf <- Vectorize(function(wL){
  BwLf1 <- Vectorize(function(w)BwLf(w, wL))
  res <- uniroot(BwLf1, c(wL, 1), tol=.Machine$double.eps)$root
  return(res)
})

# averB for invader
averBif <- function(wLi, wLr,
                    a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                    pe=-1.58*10^-3, b=4.38, h2=l*LAI/nZ*p/1000, kxmax=5, c=2.64, d=3.54,
                    gamma=1/((MAP/365/k)/1000)*nZ){
  wLLr <- wLLf(wLr)
  wLLi <- wLLf(wLi)
  
  gswLfr <- Vectorize(function(w)gswLf(w, wLr))
  gswLfi <- Vectorize(function(w)gswLf(w, wLi))
  
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

optwLif <- Vectorize(function(wLr){
  averBif1 <- Vectorize(function(wLi)averBif(wLi, wLr))
  optwLi <- optimize(averBif1, c(0.16, 0.19), tol=.Machine$double.eps, maximum=T)
  res <- optwLi$maximum-wLr
  return(res)
})

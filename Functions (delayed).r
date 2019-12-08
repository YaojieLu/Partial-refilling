
# ps(w)
psf <- function(w, pe=-1.58*10^-3, b=4.38)pe*w^(-b)

# wf(ps)
wf <- function(ps, pe=-1.58*10^-3, b=4.38)(ps/pe)^(-1/b)

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
  res <- uniroot(f1, c(-20, 0), tol=.Machine$double.eps)$root
  return(res)
})

# modified gsmax = family ESS
gsmaxfm <- function(w, wL,
                    a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=2, h=l*a*LAI/nZ*p, VPD=0.02,
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
pxminfm <- function(w, wL,
                    a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=2, h=l*a*LAI/nZ*p, VPD=0.02,
                    h2=l*LAI/nZ*p/1000, kxmax=5){
  # modified PLC
  PLCfm <- function(x)PLCf(pxL)-(PLCf(pxL)-PLCf(x))*pkx
  # modified xylem conductance function
  kxfm <- function(x)kxmax*(1-PLCfm(x))
  
  f1 <- function(x)(ps-x)*h2*kxfm(x)/(h*VPD)
  
  pxL <- psf(wL)
  ps <- psf(w)
  res <- ifelse(pxL<ps, optimize(f1, c(pxL, ps), tol=.Machine$double.eps, maximum=T)$maximum, ps)
  return(res)
}

# xylem water potential function
pxf <- function(w, gs, wL,
                a=1.6, LAI=2, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
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
Af <- function(gs, ca=400, Vcmax=50, cp=30, Km=703, Rd=1, LAI=2)LAI*1/2*(Vcmax+(Km+ca)*gs-Rd-((Vcmax)^2+2*Vcmax*(Km-ca+2*cp)*gs+((ca+Km)*gs+Rd)^2-2*Rd*Vcmax)^(1/2))

# family ESS gs(ps)
gswLpsf <- function(ps, wL){
  w <- wf(ps)
  res <- gsmaxfm(w, wL)
  return(res)
}

# family ESS A(w)
AwLf <- function(w, wL)Af(gsmaxfm(w, wL))

# ESS g1(ps)
ESSg1psf <- Vectorize(function(ps, wL, VPD=0.02, a=1.6){
  f1 <- function(w)psf(w)-ps
  w <- uniroot(f1, c(0.001, 1), tol=.Machine$double.eps)$root
  res <- sqrt(VPD*100)*(ca*gsmaxfm(w, wL)/(a*AwLf(w, wL))-1)
  return(res)
})

# averA for invader
averAif <- function(wLi, wLr,
                    a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=2, h=l*a*LAI/nZ*p, VPD=0.02,
                    pe=-1.58*10^-3, b=4.38, h2=l*LAI/nZ*p/1000, kxmax=5,
                    gamma=1/((MAP/365/k)/1000)*nZ){
  
  gswLfr <- Vectorize(function(w)ifelse(w<wLr, 0, gsmaxfm(w, wLr)))
  gswLfi <- Vectorize(function(w)ifelse(w<wLi, 0, gsmaxfm(w, wLi)))
  
  Evf <- function(w)h*VPD*gswLfr(w)
  Lf <- function(w)Evf(w)+w/100
  rLf <- function(w)1/Lf(w)
  integralrLf <- Vectorize(function(w)integrate(rLf, w, 1, rel.tol=.Machine$double.eps^0.25)$value)
  fnoc <- function(w)1/Lf(w)*exp(-gamma*w-k*integralrLf(w))
  
  f1 <- Vectorize(function(w)Af(gswLfi(w))*fnoc(w))
  res <- integrate(f1, wLi, 1, rel.tol=.Machine$double.eps^0.25)$value
  #message(wLr, " ", wLi, " ", res)
  return(res)
}

optwLif <- Vectorize(function(wLr){
  averAif1 <- Vectorize(function(wLi)averAif(wLi, wLr))
  optwLi <- optimize(averAif1, c(0.1, 0.3), tol=.Machine$double.eps^0.25, maximum=T)
  res <- optwLi$maximum-wLr
  message(wLr, " ", optwLi$maximum)
  return(res)
})

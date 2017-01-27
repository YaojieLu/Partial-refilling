
# ps(w)
psf <- function(w, pe=-1.58*10^-3, b=4.38)pe*w^(-b)

# the original PLC(px)
PLCf <- function(px, c=2.64, d=3.54)1-exp(-(-px/d)^c)

# modified kx(px)
kxfm <- function(px, wL, kxmax=5){
  # modified PLC(px)
  PLCfm <- function(px)PLCf(pxL)-(PLCf(pxL)-PLCf(px))*pkx
  
  pxL <- psf(wL)
  res <- kxmax*(1-PLCfm(px))
  return(res)
}

# modified pxmin(w)
pxminfm <- function(w, wL){
  f1 <- function(px)(ps-px)*kxfm(px, wL)
  
  ps <- psf(w)
  pxL <- psf(wL)
  res <- ifelse(pxL<ps, optimize(f1, c(pxL, ps), tol=.Machine$double.eps, maximum=T)$maximum, ps)
  return(res)
}

# xylem water potential
pxf <- function(w, gs, wL,
                a=1.6, LAI=1, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                h2=l*LAI/nZ*p/1000){
  f1 <- function(px)(ps-px)*h2*kxfm(px, wL)-h*VPD*gs
  
  ps <- psf(w)
  pxmin <- pxminfm(w, wL)
  res <- ifelse(pxmin<ps, uniroot(f1, c(pxmin, ps), tol=.Machine$double.eps)$root, ps)
  return(res)
}

# modified gsmax(w)
gsmaxfm <- function(w, wL,
                    a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                    h2=l*LAI/nZ*p/1000){
  # modified gs(px)
  f1 <- function(px)(ps-px)*h2*kxfm(px, wL)/(h*VPD)
  
  ps <- psf(w)
  pxL <- psf(wL)
  res <- ifelse(pxL<ps, optimize(f1, c(pxL, ps), tol=.Machine$double.eps, maximum=T)$objective, 0)
  return(res)
}

# Af(gs)
Af <- function(gs, ca=400, Vcmax=50, cp=30, Km=703, Rd=1, LAI=1)LAI*1/2*(Vcmax+(Km+ca)*gs-Rd-((Vcmax)^2+2*Vcmax*(Km-ca+2*cp)*gs+((ca+Km)*gs+Rd)^2-2*Rd*Vcmax)^(1/2))

# modified mf(w, gs)
mfm <- function(w, gs, wL,
                a=1.6, LAI=1, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                h2=l*LAI/nZ*p/1000, kxmax=5, c=2.64, d=3.54, h3=10){
  # modified PLC(px)
  PLCfm <- function(px)PLCf(pxL)-(PLCf(pxL)-PLCf(px))*pkx
  
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

# test family ESS
tf <- function(w, wL,
               ca=400, Vcmax=50, cp=30, Km=703, Rd=1, LAI=1,
               a=1.6, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
               h2=l*LAI/nZ*p/1000, kxmax=5, c=2.64, d=3.54, h3=10){
  
  ps <- psf(w)
  pxL <- psf(wL)
  PLCmax <- PLCf(pxL)
  # modified PLC
  PLCfm <- function(px)PLCmax-(PLCmax-PLCf(px))*pkx
  # modified xylem conductance function
  kxfm <- function(px)kxmax*(1-PLCfm(px))
  # minimum xylem water potential function for given w
  pxminf <- function(w){
    ps <- psf(w)
    f1 <- function(px)(ps-px)*h2*kxfm(px)
    res <- ifelse(pxL<ps, optimize(f1, c(pxL, ps), tol=.Machine$double.eps, maximum=T)$maximum, ps)
    return(res)
  }
  # xylem water potential function
  pxf <- function(w, gs){
    ps <- psf(w)
    pxmin <- pxminf(w)
    f1 <- function(px)((ps-px)*h2*kxfm(px)-h*VPD*gs)^2
    res <- ifelse(pxmin<ps, optimize(f1, c(pxmin, ps), tol=.Machine$double.eps)$minimum, ps)
    return(res)
  }
  
  dAdgsf <- function(gs)(1/2)*LAI*(ca+Km+((-ca^2)*gs-gs*Km^2-Km*Rd-2*cp*Vcmax-Km*Vcmax+ca*(-2*gs*Km-Rd+Vcmax))/sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))
  f1 <- function(px)(-h3*c*exp(-(-px/d)^c)*(-px/d)^(c-1)/d*pkx)*(-((exp((-(px/d))^c)*h*px*VPD)/(h2*kxmax*(exp((-(px/d))^c)*(-1+pkx)*(-1+PLCmax)*px+pkx*(px+c*ps*(-(px/d))^c-c*px*(-(px/d))^c)))))
  f2 <- function(gs)dAdgsf(gs)-f1(pxf(w, gs))
  
  res <- uniroot(f2, c(0, gsmaxfm(w, wL)), tol=.Machine$double.eps)
  return(res$root)
}

# ps(w)
psf <- function(w, pe=-1.58*10^-3, b=4.38)pe*w^(-b)

# the original PLC(px)
PLCf <- function(px, c=2.64, d=3.54)1-exp(-(-px/d)^c)

# modified gsmax
gsmaxfm <- function(w, wL,
                    a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                    h2=l*LAI/nZ*p/1000, kxmax=5, c=2.64, d=3.54){
  
  pxL <- psf(wL)
  # modified PLC
  PLCfm <- function(x)PLCf(pxL)-(PLCf(pxL)-PLCf(x))*pkx
  # modified xylem conductance function
  kxfm <- function(x)kxmax*(1-PLCfm(x))
  
  ps <- psf(w)
  f1 <- function(x)(ps-x)*h2*kxfm(x)/(h*VPD)
  
  res <- ifelse(pxL<ps, optimize(f1, c(pxL, ps), tol=.Machine$double.eps, maximum=T)$objective, 0)
  return(res)
}

# when gs reaches its max at px=pxL
ff <- function(wL,
               a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
               h2=l*LAI/nZ*p/1000, kxmax=5, c=2.64, d=3.54){
  
  f1 <- function(w){
    ps <- psf(w)
    pxL <- psf(wL)
    res <- abs(((-c)*ps*(-(pxL/d))^c+pxL*(-1+c*(-(pxL/d))^c))/(exp((-(pxL/d))^c)*pxL))
  }
  
  res <- optimize(f1, c(wL, 1), tol=.Machine$double.eps)$minimum
  return(res)
}

# Af(gs)
Af <- function(gs, ca=400, Vcmax=50, cp=30, Km=703, Rd=1, LAI=1)LAI*1/2*(Vcmax+(Km+ca)*gs-Rd-((Vcmax)^2+2*Vcmax*(Km-ca+2*cp)*gs+((ca+Km)*gs+Rd)^2-2*Rd*Vcmax)^(1/2))

# modified mf(w, gs)
mfm <- function(w, gs, wL,
                a=1.6, LAI=1, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                h2=l*LAI/nZ*p/1000, kxmax=5, c=2.64, d=3.54, h3=10){
  
  pxL <- psf(wL)
  # modified PLC
  PLCfm <- function(x)PLCf(pxL)-(PLCf(pxL)-PLCf(x))*pkx
  # modified xylem conductance function
  kxfm <- function(x)kxmax*(1-PLCfm(x))
  # minimum xylem water potential function for given w
  pxminf <- function(w){
    ps <- psf(w)
    f1 <- function(x)(ps-x)*h2*kxfm(x)
    res <- ifelse(pxL<ps, optimize(f1, c(pxL, ps), tol=.Machine$double.eps, maximum=T)$maximum, ps)
    return(res)
  }
  # xylem water potential function
  pxf <- function(w, gs){
    ps <- psf(w)
    pxmin <- pxminf(w)
    f1 <- function(x)((ps-x)*h2*kxfm(x)-h*VPD*gs)^2
    res <- ifelse(pxmin<ps, optimize(f1, c(pxmin, ps), tol=.Machine$double.eps)$minimum, ps)
    return(res)
  }
  
  px <- pxf(w, gs)
  res <- h3*(PLCfm(px)-PLCfm(0))
  return(res)
}

# modified B(w, gs)
Bfm <- function(w, gs, wL)Af(gs)-mfm(w, gs, wL)

# switch point
spf <- function(w, wL,
                ca=400, Vcmax=50, cp=30, Km=703, Rd=1, LAI=1,
                a=1.6, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                h2=l*LAI/nZ*p/1000, kxmax=5, c=2.64, d=3.54, h3=10){
  
  ps <- psf(w)
  pxL <- psf(wL)
  
  dAdgsf <- function(gs)(1/2)*LAI*(ca+Km+((-ca^2)*gs-gs*Km^2-Km*Rd-2*cp*Vcmax-Km*Vcmax+ca*(-2*gs*Km-Rd+Vcmax))/sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))
  f1 <- function(px)h3*c*exp(-(-px/d)^c)*(-px/d)^(c-1)/d*pkx*((exp((-(px/d))^c)*h*px*VPD)/(h2*kxmax*px+c*h2*kxmax*ps*(-(px/d))^c-c*h2*kxmax*px*(-(px/d))^c))
  
  #gsmax <- gsmaxfm(w, wL)
  #res <- dAdgsf(gsmax)-f1(pxL)
  gspxL <- (ps-pxL)*h2*kxmax*exp(-(-pxL/d)^c)/h/VPD
  res <- dAdgsf(gspxL)-f1(pxL)
  return(res)
}

# family ESS
gswLf <- function(w, wL){
  Bfm1 <- function(gs)Bfm(w, gs, wL)
  gsmaxfm1 <- function(w)gsmaxfm(w, wL)
  res <- ifelse(0<gsmaxfm1(w), optimize(Bfm1, c(0, gsmaxfm1(w)), tol=.Machine$double.eps, maximum=T)$maximum, 0)
  return(res)
}

# famile ESS B(w)
BwLf <- function(w, wL)Bfm(w, gswLf(w, wL), wL)

# LHS endpoint
wLLf <- Vectorize(function(wL){
  Bfm1 <- Vectorize(function(w)Bfm(w, gswLf(w, wL), wL))
  res <- uniroot(Bfm1, c(wL, 1), tol=.Machine$double.eps)$root
  return(res)
})

# averB for invader
averBif <- function(wLi, wLr,
                    a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                    pe=-1.58*10^-3, b=4.38, h2=h/1000, kxmax=5, c=2.64, d=3.54,
                    gamma=1/((MAP/365/k)/1000)*nZ){
  
  wLLr <- wLLf(wLr)
  wLLi <- wLLf(wLi)
  gswLfr <- Vectorize(function(w)ifelse(w<wLLr, 0, gswLf(w, wLr)))
  gswLfi <- Vectorize(function(w)ifelse(w<wLLi, 0, gswLf(w, wLi)))
  
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

optwLf <- Vectorize(function(wLr){
  averBif1 <- Vectorize(function(wLi)averBif(wLi, wLr))
  res <- optimize(averBif1, c(0.1, 0.3), tol=.Machine$double.eps, maximum=T)$maximum
  message(wLr)
  return(res)
})

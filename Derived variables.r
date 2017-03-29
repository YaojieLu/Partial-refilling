
options(digits=22)
source("Functions.r")
wLdata <- read.csv("Data/optwL.csv")

df <- data.frame(P50=numeric(), pxmin=numeric(), pxgs50=numeric(), slope=numeric())

for(i in 1:nrow(wLdata)){
  
  ca <- wLdata$ca[i]
  k <- wLdata$k[i]
  MAP <- wLdata$MAP[i]
  c <- wLdata$c[i]
  d <- wLdata$d[i]
  h3 <- wLdata$h3[i]
  pkx <- wLdata$pkx[i]
  wL <- wLdata$optwL[i]
  
  wLL <- wLLf(wL)
  g1 <- gswLf1(1, wL)
  f1 <- function(w)gswLf1(w, wL)-g1*0.5
  wgs50 <- uniroot(f1, c(wLL, 1), tol=.Machine$double.eps)$root
  
  df[i, 1] <- P50f(d)
  df[i, 2] <- pxf(wLL, gswLf1(wLL, wL), wL)
  df[i, 3] <- pxf(wgs50, gswLf1(wgs50, wL), wL)
  df[i, 4] <- (df[i, 2]-pxf(1, gswLf1(1, wL), wL))/(psf(wLL)-psf(1))
}

data <- cbind(wLdata, df)
write.csv(data, "Data/Derived variables.csv", row.names = FALSE)

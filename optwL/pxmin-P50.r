
options(digits=22)
source("optwL/optwL - Functions.r")

# Results
ca <- 400
k <- 0.05
MAP <- 1825
pkx <- 0.75
c <- 2.64

h3 <- 25
dSA <- seq(1, 10, by=1)
wLdata <- c(0.264470336906818, 0.218367879591088, 0.195374065471340, 0.180600369073480, 0.170034602904892, 0.161877850407940, 0.155301800064056, 0.149859211988843, 0.145244256511873, 0.141216842521499)
df <- data.frame(P50=rep(0, length(dSA)), pxmin=rep(0, length(dSA)))

for(i in 1:length(dSA)){
  d <- dSA[i]
  wL <- wLdata[i]
  wLL <- wLLf(wL)
  df[i, 1] <- P50f(d)
  df[i, 2] <- pxf(wLL, gswLf1(wLL, wL), wL)
}
df

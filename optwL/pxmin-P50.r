
options(digits=20)
source("optwL/optwL - Functions.r")

# Results
ca <- 400
k <- 0.05
MAP <- 1825
pkx <- 0.5
c <- 2.64

wLdata <- c(0.206910328, 0.184815566, 0.170278237, 0.159746981, 0.151582936, 0.144992154, 0.139475231, 0.134854655, 0.130846915)
ddata <- seq(2, 10, by=1)

df <- data.frame(pxmin=rep(0, length(ddata)), P50=rep(0, length(ddata)))

for(i in 1:length(ddata)){
  d <- ddata[i]
  wL <- wLdata[i]
  df[i, 1] <- pxminfm(wL, wL)
  df[i, 2] <- P50f(wL)
}
plot(df[, 2], df[, 1], xlim=c(-15, 0), ylim=c(-15, 0))

# No refilling
dfno <- data.frame(pxmin=c(-1.840809841, -2.921439074, -4.076371045, -5.288027902, -6.544306864, -7.836769373, -9.159363146, -10.50761004, -11.87809899), P50=c(-1.740750825, -2.611126238, -3.48150165, -4.351877063, -5.222252475, -6.092627888, -6.9630033, -7.833378713, -8.703754125))
points(dfno[, 2], dfno[, 1], col="red")

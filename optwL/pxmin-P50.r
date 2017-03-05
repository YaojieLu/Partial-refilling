
options(digits=20)
source("optwL/optwL - Functions.r")

# Results
ca <- 400
k <- 0.05
MAP <- 1825
pkx <- 0.5
c <- 2.64

h3 <- 10
dSA <- seq(2, 10, by=1)
wLdata <- c(0.206910328306461, 0.184815566462289, 0.170278236900065, 0.159746981368831, 0.151582935969734, 0.144992153977275, 0.139475230518173, 0.134854655191347, 0.130846914876583)
df <- data.frame(pxmin=rep(0, length(dSA)), P50=rep(0, length(dSA)))

for(i in 1:length(dSA)){
  d <- dSA[i]
  wL <- wLdata[i]
  df[i, 1] <- pxminfm(wL, wL)
  df[i, 2] <- P50f(wL)
}
df

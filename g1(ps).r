
source("Functions.r")
data <- read.csv("Data/Derived variables.csv")

# Initializing
ca <- 400
c <- 2.64
d <- 5

h3SA <- c(25, 100, 25, 25)
pkxSA <- c(0.5, 0.5, 0.25, 0.75)
wLSA <- c(as.numeric(subset(data, k==0.05 & MAP==1825 & d==5 & h3==h3SA[1] & pkx==pkxSA[1], select=c("optwL"))),
          as.numeric(subset(data, k==0.05 & MAP==1825 & d==5 & h3==h3SA[2] & pkx==pkxSA[2], select=c("optwL"))),
          as.numeric(subset(data, k==0.05 & MAP==1825 & d==5 & h3==h3SA[3] & pkx==pkxSA[3], select=c("optwL"))),
          as.numeric(subset(data, k==0.05 & MAP==1825 & d==5 & h3==h3SA[4] & pkx==pkxSA[4], select=c("optwL"))))
res1 <- data.frame(h3=c(rep(h3SA[1], 101), rep(h3SA[2], 101), rep(h3SA[3], 101), rep(h3SA[4], 101)),
                   pkx=c(rep(pkxSA[1], 101), rep(pkxSA[2], 101), rep(pkxSA[3], 101), rep(pkxSA[4], 101)))
res2 <- data.frame(ps=numeric(), g1=numeric())

for(i in 1:length(wLSA)){
  h3 <- h3SA[i]
  pkx <- pkxSA[i]
  wL <- wLSA[i]
  wLL <- wLLf(wL)
  pxL <- psf(wLL)
  f1 <- function(ps)ESSg1psf(ps, wL)
  x <- seq(pxL, psf(1), by=(psf(1)-pxL)/100)
  y <- f1(x)
  datat <- data.frame(x, y)
  res2 <- rbind(res2, datat)
}
res <- cbind(res1, res2)
colnames(res) <- c("h3", "pkx", "ps", "g1")
write.csv(res, "Data/g1(ps).csv", row.names = FALSE)

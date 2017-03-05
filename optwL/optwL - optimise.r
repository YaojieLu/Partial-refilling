
options(digits=20)
source("optwL/optwL - Functions.r")
data <- read.csv("data/optwL - pkx.csv")

# Results
ca <- 400
k <- 0.05
MAP <- 1825
pkx <- 0.5
c <- 2.64

h3SA <- 1
h3 <- h3SA
df <- data.frame(wL0=subset(data, h3==as.numeric(h3SA), select=c("No.refilling")), wL100=subset(data, h3==h3SA, select=c("Perfect.refilling")))
optwL <- vector("list")

for(i in 1:10){
  d <- i
  optwL[[i]] <- try(uniroot(optwLif, sort(df[i, ]), tol=.Machine$double.eps^0.25))
}
optwL

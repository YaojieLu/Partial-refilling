
options(digits=20)
source("Functions.r")
wLdata <- read.csv("Data/optwL.csv")

ca <- 400
c <- 2.64

df <- data.frame(averA=numeric(), averE=numeric())
for(i in 1:1){
  k <- wLdata[i, 2]
  MAP <- wLdata[i, 3]
  h3 <- wLdata[i, 5]
  d <- wLdata[i, 6]
  pkx <- wLdata[i, 7]
  wL <- wLdata[i, 8]
  df <- averBf(wL)
}

#data <- cbind(wLdata, df)
#write.csv(data, "Data/Derived variables.csv", row.names = FALSE)

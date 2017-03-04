
options(digits=20)
source("optwL/optwL - Functions.r")

# Results
ca <- 400
k <- 0.05
MAP <- 1825
pkx <- 0.5

c <- 2.64
d <- 10

SAr <- c(0.175, 0.177, 0.176)
SAi <- c(0.173, 0.176, 0.179)
df <- as.vector(expand.grid(wLi=SAi, wLr=SAr))
df[, "averBi"] <- numeric(length=nrow(df))

for(i in 1:nrow(df)){
  df[i, 3] <- averBif(wLi=df[i, 1], wLr=df[i, 2])
}


options(digits=20)
source("Functions.r")

# Results
ca <- 400
k <- 0.05
MAP <- 1825
pkx <- 0.5

SAr <- seq(0.1776, 0.178, by=0.0001)
SAi <- seq(0.17, 0.19, by=0.003)
df <- as.vector(expand.grid(wLi=SAi, wLr=SAr))
df[, "averBi"] <- numeric(length=nrow(df))

for(i in 1:nrow(df)){
  df[i, 3] <- averBif(wLi=df[i, 1], wLr=df[i, 2])
}

write.csv(df, "Results/Invasion analysis by=0.0001,0.003.csv", row.names = FALSE)

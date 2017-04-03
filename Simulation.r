
options(digits=20)
source("Functions.r")

# Results
ca <- 400
c <- 2.64

k <- 0.05
MAP <- 1825

pkx <- 0.5
h3 <- 25
d <- 5

SA <- 
  df <- as.vector(expand.grid(wLi=seq(0.16, 0.18, by=0.001), wLr=seq(0.16, 0.18, by=0.001)))
df[, "averBi"] <- numeric(length=nrow(df))

for(i in 1:nrow(df)){
  df[i, 3] <- averBif(wLi=df[i, 1], wLr=df[i, 2])
}

#write.csv(df, "Data/Invasion analysis.csv", row.names = FALSE)

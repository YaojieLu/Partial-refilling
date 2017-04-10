
options(digits=20)
source("tf.r")

# Results
ca <- 400
c <- 2.64

k <- 0.05
MAP <- 1825
h3 <- 1
pkx <- 0.5
d <- 1

wL <- 0.24
x <- seq(0.2, 0.3, by=0.01)
y <- numeric()
for(i in 1:length(x)){
  y[i] <- averBif(wLi=x[i], wLr=wL)
}
plot(x, y)
abline(v=wL)

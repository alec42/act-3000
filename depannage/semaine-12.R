# Code R presente lors du 24 novembre 2017
# Christopher Blier-Wong
# ACT-3000 Theorie du risque

# Plan du depannage -------------------------------------------------------

# Méthode des rectangles
# Calcul du rho de Spearman empirique

# Méthode des rectangles --------------------------------------------------

# On compare premièrement les graphiques
DessinRectangle <- function(m){
  par(mfrow = c(1, 2))
  
  s <- 5
  pas <- s / 2 ** m
  x.lower <- seq(0, s - pas, by = pas)
  x.upper <- seq(0, s, by = pas)
  
  plot(0:10, 0:10, main = paste("Méthode Lower m =", m), type = "n", xaxs = "i", yaxs = "i", 
       xlab = expression(X[1]), ylab = expression(X[2]))
  axis(1, s, expression(italic(s)))
  axis(2, s, expression(italic(s)))
  
  for(i in x.lower){
    polygon(c(i, i, i + pas, i + pas), c(0, s - i - pas, s - i - pas, 0), col = "red")
  }
  
  segments(0, s, s, 0, col = "blue")
  
  plot(0:10, 0:10, main = paste("Méthode Upper m =", m), type = "n", xaxs = "i", yaxs = "i", 
       xlab = expression(X[1]), ylab = expression(X[2]))
  axis(1, 5, expression(italic(s)))
  axis(2, 5, expression(italic(s)))
  
  for(i in x.upper){
    polygon(c(i, i, i + pas, i + pas), c(0, s - i  ,s - i, 0), col = "red")
  }
  segments(0, s, s, 0, col = "blue")
  
  par(mfrow=c(1, 1))
}

DessinRectangle(1)
DessinRectangle(2)
DessinRectangle(3)
DessinRectangle(4)
DessinRectangle(5)
DessinRectangle(6)
DessinRectangle(10)

# Retour sur la question 10.1.1 # 4 

RepartitionCopule <- function(u1, u2) {
  u1 + u2 - 1 + ((1 - u1) ** (- 2) + (1 - u2) ** (- 2) - 1) ** (- 1 / 2)
}

Fy1 <- function(y1) {
  1 - (20 / (20 + y1)) ** 3
}

Fy2 <- function(y2) {
  pgamma(y2, shape = 2, rate = 1 / 10)
}

Fy1y2 <- function(y1, y2) {
  RepartitionCopule(Fy1(y1), Fy2(y2))
}

CalcuerRepartition <- function(m, s){
  pas <- s / (2 ** m)  
  ylower <- seq(pas, s - pas, by = pas)
  yupper <- seq(pas, s, by = pas)
  a.lower <- sum(sapply(ylower, function(t) Fy1y2(t, s - t) - Fy1y2(t - pas, s - t)))
  a.upper <- sum(sapply(yupper, function(t) Fy1y2(t, s - t + pas) - Fy1y2(t - pas, s - t  + pas)))
  print(paste0("Lower: ", a.lower, ". Upper: ",  a.upper, "."))
}

CalcuerRepartition(1, 60)
CalcuerRepartition(2, 60)
CalcuerRepartition(3, 60)
CalcuerRepartition(4, 60)
CalcuerRepartition(5, 60)
CalcuerRepartition(6, 60)
CalcuerRepartition(10, 60)


# Rho de Spearman ---------------------------------------------------------

# Marginales log-normale et differentes copules

nsim <- 1000
alphabeta <- c(.3, .3, 1-.3-.3)
vSup <- matrix(runif(nsim), nsim, 2)
vInf <- matrix(runif(nsim), nsim, 2)
vInf[,2]<-1 - vInf[,1]
vInd <- matrix(runif(2 * nsim), nsim, 2, byrow = T)
fb <- sample(c(1, 2, 3), 1000, replace = TRUE, prob = alphabeta)
vV <- matrix(numeric(), nsim, 2)

for(i in 1:nsim) {
  if (fb[i] == 1) {
    vV[i, ] <- vInf[i, ]
  } else if (fb[i] == 2) {
    vV[i, ] <- vSup[i,]
  } else {
    vV[i, ] <- vInd[i, ]
  }
}

vX <- cbind(qlnorm(vV[, 1], log(10), 0.2), qlnorm(vV[, 2], log(20), 0.3))

(mean(rank(vX[, 1]) * rank(vX[, 2])) - mean(rank(vX[, 1]) * mean(rank(vX[, 2])))) / 
  sqrt(var(rank(vX[, 1])) * var(rank(vX[, 2])))

cor(vX, method = "spearman")[2]

diff(alphabeta)[1]

# Copule EFGM

nsim <- 10000
alpha <- 0.5
vV <- matrix(runif(nsim * 2), nsim, 2, byrow = T)
W1 <- alpha * (2 * vV[ ,1] - 1) - 1
W2 <- (1 - alpha * (2 * vV[ ,1] - 1)) ** 2 + 4 * alpha * vV[ ,2] * (2 * vV[ ,1] - 1)
vU <- cbind(vV[ ,1], 2 * vV[ ,2] / (sqrt(W2) - W1))

vX <- cbind(qlnorm(vU[, 1], log(10), 0.9), qlnorm(vU[, 2], log(200), 0.3))

(mean(rank(vX[, 1]) * rank(vX[, 2])) - mean(rank(vX[, 1]) * mean(rank(vX[, 2])))) / 
  sqrt(var(rank(vX[, 1])) * var(rank(vX[, 2])))

cor(vX, method = "spearman")[2]

alpha / 3

# Copule de Clayton

alph <- 10
vV <- matrix(runif(nsim * 3), nsim, 3, byrow = T)
vTheta <- qgamma(vV[, 1], 1 / alph, 1)
vY <- sapply(1:2, function(t) qexp(vV[, t + 1], vTheta))
vU <- (1 + vY) ** (-1 / alph)

vX <- cbind(qlnorm(vU[, 1], log(1), 0.1), qlnorm(vU[, 2], log(2000), 0.9))

(mean(rank(vX[, 1]) * rank(vX[, 2])) - mean(rank(vX[, 1]) * mean(rank(vX[, 2])))) / 
  sqrt(var(rank(vX[, 1])) * var(rank(vX[, 2])))

cor(vX, method = "spearman")[2]

# Estimation avec la copule EFGM ------------------------------------------

set.seed(20171123)
nsim <- 2000
alpha <- - 0.5
vV <- matrix(runif(nsim * 2), nsim, 2, byrow = T)
W1 <- alpha * (2 * vV[ ,1] - 1) - 1
W2 <- (1 - alpha * (2 * vV[ ,1] - 1)) ** 2 + 4 * alpha * vV[ ,2] * (2 * vV[ ,1] - 1)
vU <- cbind(vV[ ,1], 2 * vV[ ,2] / (sqrt(W2) - W1))
plot(vU, xlab = expression(U[1]), ylab = expression(U[2])) 

# Estimation du paramètre de la copule

# Méthode des moments: 

# Rho de Spearman
# On sait que rho = alpha / 3
# Alors alpha = 3 * rho de Spearman

3 * cor(vU, method = "spearman")[2]

# Tau de Kendall
# On sait que tau = 2 * alpha / 9
# Alors alpha = 9 * tau / 2

9 * cor(vU, method = "kendall")[2] / 2

# Code R presente lors du 24 novembre 2017
# Christopher Blier-Wong
# ACT-3000 Theorie du risque

# Plan du depannage -------------------------------------------------------

# Méthode des rectangles


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


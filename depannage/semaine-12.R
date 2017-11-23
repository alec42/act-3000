# Code R presente lors du 24 novembre 2017
# Christopher Blier-Wong
# ACT-3000 Theorie du risque

# Plan du depannage -------------------------------------------------------

# Méthode des rectangles
# Estimation 



# Méthode des rectangles --------------------------------------------------

# On compare premièrement les graphiques
DessinRectangle <- function(m){
  par(mfrow = c(1, 2))
  
  s <- 5
  pas <- s / 2 ** m
  x.lower <- seq(0, s - pas, by = pas)
  x.upper <- seq(0, s, by = pas)
  
  plot(0:10, 0:10, main = "Méthode Lower", type = "n", xaxs = "i", yaxs = "i", 
       xlab = expression(X[1]), ylab = expression(X[2]))
  axis(1, s, expression(italic(s)))
  axis(2, s, expression(italic(s)))
  
  for(i in x.lower){
    polygon(c(i, i, i + pas, i + pas), c(0, s - i - pas, s - i - pas, 0), col = "red")
  }
  
  segments(0, s, s, 0, col = "blue")
  
  plot(0:10, 0:10, main = "Méthode Upper", type = "n", xaxs = "i", yaxs = "i", 
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

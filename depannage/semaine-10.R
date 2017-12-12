# Code R presente lors du 10 novembre 2017
# Christopher Blier-Wong
# ACT-3000 Theorie du risque

# Plan du depannage -------------------------------------------------------

# Exemples avec les methodes de discretisation
# Inversion de lois bivariees 
# Creation de lois bivariees avec les copules
# Exemples de simulation de copules

# Exemples discretisation -------------------------------------------------

# La methode de discretisation upper est simplement une translation 
# de un de la methode lower. 

# Voici un exemple avec la distribution gamma(2, 10)
# On discretise la fonction de repartition par pas de 
# h entre x.min et x.max

h <- 0.1
x.min <- 0
x.max <- 1000

repartition.discrete <- pgamma(seq(from = x.min, to = x.max + h, by = h), 
                               shape = 2, rate = 1 / 10)

# Si h = 1 x.min = 0, on peut simplement faire 
#repartition.discrete <- pgamma(0:(x.max + 1), shape = 2, rate = 1 / 10)

densite.upper <- diff(repartition.discrete)
densite.lower <- c(0, densite.upper[- x.max / h])

# Comparaison d'esperances

sum(seq(x.min, x.max, h) * densite.upper)  # Upper sous-estime
sum(seq(x.min, x.max, h) * densite.lower)  # Lower sur-estime
2 * 10  # Theorique
2 * 10 + c(- h / 2, h / 2)  # Les valeurs sont sous ou sur estimes 
                            # par h / 2

# Exemple avec fft

h <- 0.1
a <- 2 ** 14
phi.b <- fft(c(densite.upper, numeric(a - x.max / h)))
phi.s <- exp(10 * (phi.b - 1))
fs.upper <- Re(fft(phi.s, inverse = TRUE)) / a

phi.b <- fft(c(densite.lower, numeric(a - x.max / h)))
phi.s <- exp(10 * (phi.b - 1))
fs.lower <- Re(fft(phi.s, inverse = TRUE)) / a

sum(seq(0, by = h, length.out = (a + 1)) * fs.upper)
sum(seq(0, by = h, length.out = (a + 1)) * fs.lower)

# Examen paritiel info A16 numero 3

h <- 1
x.max <- 500

# (a)
fb <- c(0, diff(plnorm(0:(x.max + h), log(10) - 0.36 / 2, 0.6)))
fb[c(1, 11)]  # Reponse
fb[16]  # Verification

# (b)

panjer.nbinom1<-function(rr, qq, ff, smax)
{
  # Algorithme de Panjer
  # Cas Binomiale negative 1
  # Loi discrete pour B
  aa<-1-qq
  bb<-aa*(rr-1)
  ll<-length(ff)
  ffs<-(qq/(1-(1-qq)*ff[1]))^rr
  ff<-c(ff,rep(0,smax-ll+1))
  for (i in 1 :smax)
  {
    j<-i+1
    ffs<-c(ffs,(1/(1-aa*ff[1]))*sum(ff[2 :j]*ffs[i :1]*(bb*(1 :i)/i+aa)))
  }
  return(ffs)
}


fx <- panjer.nbinom1(1.5, 1 / 3, fb, 600)
fx[c(1, 21)]  # Reponse
fx[16]  # Verification

Fx <- cumsum(fx)
Fx[51]  # Reponse
Fx[61]  # Verification

stop.loss <- sum(pmax(0:(length(fx) - 1) - 50, 0) * fx)

50 + 1 / (1 - Fx[51]) * stop.loss

# Exercice copule 10.1.1 # 4 ----------------------------------------------

DensiteCopule <- function(u1, u2) {
  u1 + u2 - 1 + ((1 - u1) ** (- 2) + (1 - u2) ** (- 2) - 1) ** (- 1 / 2)
}

Fy1 <- function(y1) {
  1 - (20 / (20 + y1)) ** 3
}

Fy2 <- function(y2) {
  pgamma(y2, shape = 2, rate = 1 / 10)
}

Fy1y2 <- function(y1, y2) {
  DensiteCopule(Fy1(y1), Fy2(y2))
}

# (i)

Fy1y2(15, 11) - Fy1y2(5, 11) - Fy1y2(15, 6) + Fy1y2(5, 6)

# (ii)

(1 - Fy1(50) - Fy2(30) + Fy1y2(50, 30)) / (1 - Fy2(30))

# (iii)

s <- 60
m <- 15

pas <- s / (2 ** m)

# Methode lower
ylower <- seq(pas, s - pas, by = pas)
sum(sapply(ylower, function(t) Fy1y2(t, s - t) - Fy1y2(t - pas, s - t)))

# Methode upper
yupper <- seq(pas, s, by = pas)
sum(sapply(yupper, function(t) Fy1y2(t, s - t + pas) - Fy1y2(t - pas, s - t  + pas)))

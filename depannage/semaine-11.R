# Code R presente lors du 16 novembre 2017
# Christopher Blier-Wong
# ACT-3000 Theorie du risque

# Plan du depannage -------------------------------------------------------

# Mesures de risque et copules
# Exercices traditionnels
# Anciennes questions d'examen

# Mesures de risque et copules --------------------------------------------

# Cette section montre les solutions des exercices de la section 10.5 du 
# document d'exercice, déposé sur le site du cours le 16 novembre 2017. 

library(pracma)

# 10.5.1

sapply(c(10**-3, 10**-6), function(kap) ((2 * kap ** -5 - 1) ** -(1 / 5)) / kap)
curve(((2 * x ** -5 - 1) ** -(1 / 5)) / x)
sapply(c(10 ** -3, 10 ** -6), function(kap) (2 * kap - 1 + (2 * (1 - kap) ** -5 - 1) ** -(1 / 5)) / kap)

sapply(c(1 - 10 ** -3, 1 - 10 ** -6), function(kap) (2 * kap - 1 + ((2 * (1 - kap) ** -5 - 1) ** -(1 / 5))) / kap)
# Dépendance très positive aux extrèmes.

sapply(c(1 - 10 ** -3, 1 - 10 ** -3), function(kap) ((2 * (1 - kap) ** -5 - 1) ** -(1 / 5)) / (1 - kap))

# (c) ii

alpha <- 5

C.clayton <- function(u, v) (u ** (-alpha) + v ** (-alpha) - 1) ** (-1 / alpha)
C.clayton.bar <- function(u, v) u + v - 1 + C.clayton(1 - u, 1 - v)

# U_1, U_2

Fu1u2 <- function(u, v) C.clayton(u, v)
Fu1u2.bar <- function(u, v) 1 - u - v + Fu1u2(u, v)
Eu1u2 <- quad2d(Fu1u2.bar, 0, 1, 0, 1)
covu1u2 <- Eu1u2 - 0.5 ** 2
covu1u2 / sqrt(1 / 12 * 1 / 12)

# U_1', U_2'

Fu1u2 <- function(u, v) C.clayton.bar(u, v)
Fu1u2.bar <- function(u, v) 1 - u - v + Fu1u2(u, v)
Eu1u2 <- quad2d(Fu1u2.bar, 0, 1, 0, 1)
covu1u2 <- Eu1u2 - 0.5 ** 2
covu1u2 / sqrt(1 / 12 * 1 / 12)

# (c) iii

EX <- 1
VX <- 4

F1 <- function(x) pgamma(x, 1/4, 1/4)
F2 <- function(x) pgamma(x, 1/4, 1/4)

Fx1x2 <- function(x, y) C.clayton(F1(x), F2(y))
Fx1x2.bar <- function(x, y) 1 - F1(x) - F2(y) + Fx1x2(x, y)
Ex1x2 <- quad2d(Fx1x2.bar, 0, 100, 0, 100)
covx1x2 <- Ex1x2 - EX ** 2
covx1x2 / sqrt(VX ** 2)

Fx1x2 <- function(x, y) C.clayton.bar(F1(x), F2(y))
Fx1x2.bar <- function(x, y) 1 - F1(x) - F2(y) + Fx1x2(x, y)
Ex1x2 <- quad2d(Fx1x2.bar, 0, 100, 0, 100)
covx1x2 <- Ex1x2 - EX ** 2
covx1x2 / sqrt(VX ** 2)

# (c) iv

nsim <- 1000000

vV <- matrix(runif(nsim * 3), nsim, 3, byrow = T)
vTheta <- qgamma(vV[, 1], 1 / alpha, 1)
vY <- sapply(1:2, function(t) qexp(vV[, t + 1], vTheta))
vU <- (1 + vY) ** (-1 / alpha)

(mean(vU[, 1] * vU[, 2]) - prod(colMeans(vU))) / sqrt( var(vU[, 1]) * var(vU[, 2]) )

# (c) v

vU.prime <- 1 - vU
(mean(vU.prime[, 1] * vU.prime[, 2]) - prod(colMeans(vU.prime))) / sqrt( var(vU.prime[, 1]) * var(vU.prime[, 2]) )

# (c) vi

X1X2 <- qgamma(vU, 1 / 4, 1 / 4)
X1X2.prime <- qgamma(vU.prime, 1 / 4, 1 / 4)

# (c) vii

S <- rowSums(X1X2)
S.prime <- rowSums(X1X2.prime)

# (c) viii

(VaRS <- sapply(c(0.9, 0.99, 0.999, 0.9999), function(t) sort(S)[t * nsim]))
(VaRS.prime <- sapply(c(0.9, 0.99, 0.999, 0.9999), function(t) sort(S.prime)[t * nsim]))

(TVaRS <- sapply(VaRS, function(t) mean(S[S > t])))
(TVaRS <- sapply(VaRS.prime, function(t) mean(S.prime[S.prime > t])))

# Effectuer les calculs suivants pour les deux hypothèses et pour alpha tel que rhoP = 0.5

F1 <- function(x) pgamma(x, 1/4, 1/4)
F2 <- function(x) pgamma(x, 1/4, 1/4)

rhox1x2 <- function(alpha){
  C.clayton <- function(u, v) (u ** (-alpha) + v ** (-alpha) - 1) ** (-1 / alpha)
  Fx1x2 <- function(x, y) C.clayton(F1(x), F2(y))
  Fx1x2.bar <- function(x, y) 1 - F1(x) - F2(y) + Fx1x2(x, y)
  Ex1x2 <- quad2d(Fx1x2.bar, 0, 100, 0, 100)
  covx1x2 <- Ex1x2 - EX ** 2
  covx1x2 / sqrt(VX ** 2)
}

TrouverAlpha <- function(rho){
  optimize(function(x) abs(rhox1x2(x) - rho), c(0, 10))$minimum
}

TrouverAlpha(0.5)

rhox1x2 <- function(alpha){
  C.clayton <- function(u, v) (u ** (-alpha) + v ** (-alpha) - 1) ** (-1 / alpha)
  C.clayton.bar <- function(u, v) u + v - 1 + C.clayton(1 - u, 1 - v)
  Fx1x2 <- function(x, y) C.clayton.bar(F1(x), F2(y))
  Fx1x2.bar <- function(x, y) 1 - F1(x) - F2(y) + Fx1x2(x, y)
  Ex1x2 <- quad2d(Fx1x2.bar, 0, 100, 0, 100)
  covx1x2 <- Ex1x2 - EX ** 2
  covx1x2 / sqrt(VX ** 2)
}

TrouverAlpha <- function(rho){
  optimize(function(x) abs(rhox1x2(x) - rho), c(0, 10))$minimum
}

TrouverAlpha(0.5)







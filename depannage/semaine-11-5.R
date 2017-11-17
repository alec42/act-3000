library(pracma)

# 10.5.1

# (c) ii

rstable <- function(n, phi){ 
  # algo de Marceau (2013) page 351 
  R <- runif(n, min = - pi / 2, max = pi / 2) 
  V <- rexp(n, 1) 
  A <- sin(phi * (R + pi / 2)) / cos(R) ** (1 / phi) 
  B <- (cos(phi * pi / 2 + (phi - 1) * R) / V) ** ((1 - phi) / phi) 
  A * B 
} 

alpha <- 5

C.gumbel <- function(u, v) exp(-((-log(u) ** alpha) + (-log(v)) ** alpha) ** (1 / alpha))
C.gumbel.bar <- function(u, v) u + v - 1 + C.gumbel(1 - u, 1 - v)

# U_1, U_2

Fu1u2 <- function(u, v) C.gumbel(u, v)
Fu1u2.bar <- function(u, v) 1 - u - v + Fu1u2(u, v)
Eu1u2 <- quad2d(Fu1u2.bar, 0, 1, 0, 1)
covu1u2 <- Eu1u2 - 0.5 ** 2
covu1u2 / sqrt(1 / 12 * 1 / 12)

# U_1', U_2'

Fu1u2 <- function(u, v) C.gumbel.bar(u, v)
Fu1u2.bar <- function(u, v) 1 - u - v + Fu1u2(u, v)
Eu1u2 <- quad2d(Fu1u2.bar, 0, 1, 0, 1)
covu1u2 <- Eu1u2 - 0.5 ** 2
covu1u2 / sqrt(1 / 12 * 1 / 12)

# (c) iii

EX <- 1
VX <- 4

F1 <- function(x) pgamma(x, 1/4, 1/4)
F2 <- function(x) pgamma(x, 1/4, 1/4)

Fx1x2 <- function(x, y) C.gumbel(F1(x), F2(y))
Fx1x2.bar <- function(x, y) 1 - F1(x) - F2(y) + Fx1x2(x, y)
Ex1x2 <- quad2d(Fx1x2.bar, 0, 100, 0, 100)
covx1x2 <- Ex1x2 - EX ** 2
covx1x2 / sqrt(VX ** 2)

Fx1x2 <- function(x, y) C.gumbel.bar(F1(x), F2(y))
Fx1x2.bar <- function(x, y) 1 - F1(x) - F2(y) + Fx1x2(x, y)
Ex1x2 <- quad2d(Fx1x2.bar, 0, 100, 0, 100)
covx1x2 <- Ex1x2 - EX ** 2
covx1x2 / sqrt(VX ** 2)

# (c) iv

nsim <- 1000000
set.seed(2017)
vThet <- rstable(nsim, phi = alpha ** -1) 
vV<-matrix(runif(nsim * 2), nsim, 2, byrow = TRUE)
vZ <- qexp(vV, 1)
vU <- exp(- (vZ / vThet) ** (1 / alpha)) 

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

# plot(ecdf(S))
# plot(ecdf(S.prime), add = TRUE)

# (c) viii

(VaRS <- sapply(c(0.9, 0.99, 0.999, 0.9999), function(t) sort(S)[t * nsim]))
(VaRS.prime <- sapply(c(0.9, 0.99, 0.999, 0.9999), function(t) sort(S.prime)[t * nsim]))

(TVaRS <- sapply(VaRS, function(t) mean(S[S > t])))
(TVaRS <- sapply(VaRS.prime, function(t) mean(S.prime[S.prime > t])))

# Effectuer les calculs suivants pour les deux hypothèses et pour alpha tel que rhoP = 0.5

F1 <- function(x) pgamma(x, 1/4, 1/4)
F2 <- function(x) pgamma(x, 1/4, 1/4)

rhox1x2 <- function(alpha){
  C.gumbel <- function(u, v) exp(-((-log(u) ** alpha) + (-log(v)) ** alpha) ** (1 / alpha))
  Fx1x2 <- function(x, y) C.gumbel(F1(x), F2(y))
  Fx1x2.bar <- function(x, y) 1 - F1(x) - F2(y) + Fx1x2(x, y)
  Ex1x2 <- quad2d(Fx1x2.bar, 0, 50, 0, 50)
  covx1x2 <- Ex1x2 - EX ** 2
  covx1x2 / sqrt(VX ** 2)
}

TrouverAlpha <- function(rho){
  optimize(function(x) abs(rhox1x2(x) - rho), c(0, 10))$minimum
}

(alpha <- TrouverAlpha(0.5))

sapply((1 : 10) , rhox1x2)

# (c) iv

nsim <- 1000000
set.seed(2017)
vThet <- rstable(nsim, phi = alpha ** -1) 
vV<-matrix(runif(nsim * 2), nsim, 2, byrow = TRUE)
vZ <- qexp(vV, 1)
vU <- exp(- (vZ / vThet) ** (1 / alpha)) 



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

# plot(ecdf(S))
# plot(ecdf(S.prime), add = TRUE)

# (c) viii

(VaRS <- sapply(c(0.9, 0.99, 0.999, 0.9999), function(t) sort(S)[t * nsim]))
(VaRS.prime <- sapply(c(0.9, 0.99, 0.999, 0.9999), function(t) sort(S.prime)[t * nsim]))

(TVaRS <- sapply(VaRS, function(t) mean(S[S > t])))
(TVaRS <- sapply(VaRS.prime, function(t) mean(S.prime[S.prime > t])))


rhox1x2 <- function(alpha){
  C.gumbel <- function(u, v) exp(-((-log(u) ** alpha) + (-log(v)) ** alpha) ** (1 / alpha))
  C.gumbel.bar <- function(u, v) u + v - 1 + C.gumbel(1 - u, 1 - v)
  Fx1x2 <- function(x, y) C.gumbel.bar(F1(x), F2(y))
  Fx1x2.bar <- function(x, y) 1 - F1(x) - F2(y) + Fx1x2(x, y)
  Ex1x2 <- quad2d(Fx1x2.bar, 0, 100, 0, 100)
  covx1x2 <- Ex1x2 - EX ** 2
  covx1x2 / sqrt(VX ** 2)
}

TrouverAlpha <- function(rho){
  optimize(function(x) abs(rhox1x2(x) - rho), c(0, 10))$minimum
}

sapply(1:10, rhox1x2)


# (c) iv

nsim <- 1000000
set.seed(2017)
vThet <- rstable(nsim, phi = alpha ** -1) 
vV<-matrix(runif(nsim * 2), nsim, 2, byrow = TRUE)
vZ <- qexp(vV, 1)
vU <- exp(- (vZ / vThet) ** (1 / alpha)) 

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

# plot(ecdf(S))
# plot(ecdf(S.prime), add = TRUE)

# (c) viii

(VaRS <- sapply(c(0.9, 0.99, 0.999, 0.9999), function(t) sort(S)[t * nsim]))
(VaRS.prime <- sapply(c(0.9, 0.99, 0.999, 0.9999), function(t) sort(S.prime)[t * nsim]))

(TVaRS <- sapply(VaRS, function(t) mean(S[S > t])))
(TVaRS <- sapply(VaRS.prime, function(t) mean(S.prime[S.prime > t])))


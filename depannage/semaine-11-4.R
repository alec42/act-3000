library(pracma)

# 10.5.1

# (c) ii

alpha <- 5

C.frank <- function(u, v) -1 / alpha * log(1 + (exp(-alpha * u) - 1) * (exp(-alpha * v) - 1) / (exp(-alpha) - 1))
C.frank.bar <- function(u, v) u + v - 1 + C.frank(1 - u, 1 - v)

# U_1, U_2

Fu1u2 <- function(u, v) C.frank(u, v)
Fu1u2.bar <- function(u, v) 1 - u - v + Fu1u2(u, v)
Eu1u2 <- quad2d(Fu1u2.bar, 0, 1, 0, 1)
covu1u2 <- Eu1u2 - 0.5 ** 2
covu1u2 / sqrt(1 / 12 * 1 / 12)

# U_1', U_2'

Fu1u2 <- function(u, v) C.frank.bar(u, v)
Fu1u2.bar <- function(u, v) 1 - u - v + Fu1u2(u, v)
Eu1u2 <- quad2d(Fu1u2.bar, 0, 1, 0, 1)
covu1u2 <- Eu1u2 - 0.5 ** 2
covu1u2 / sqrt(1 / 12 * 1 / 12)

# (c) iii

EX <- 1
VX <- 4

F1 <- function(x) plnorm(x, - 0.5 * log(5), sqrt(log(5)))
F2 <- function(x) plnorm(x, - 0.5 * log(5), sqrt(log(5)))


Fx1x2 <- function(x, y) C.frank(F1(x), F2(y))
Fx1x2.bar <- function(x, y) 1 - F1(x) - F2(y) + Fx1x2(x, y)
Ex1x2 <- quad2d(Fx1x2.bar, 0, 100, 0, 100)
covx1x2 <- Ex1x2 - EX ** 2
covx1x2 / sqrt(VX ** 2)

Fx1x2 <- function(x, y) C.frank.bar(F1(x), F2(y))
Fx1x2.bar <- function(x, y) 1 - F1(x) - F2(y) + Fx1x2(x, y)
Ex1x2 <- quad2d(Fx1x2.bar, 0, 100, 0, 100)
covx1x2 <- Ex1x2 - EX ** 2
covx1x2 / sqrt(VX ** 2)

# (c) iv

nsim <- 1000000
set.seed(2017)
vV <- matrix(runif(nsim * 2), nsim, 2, byrow = TRUE)
vU <- cbind(vV[ ,1], -1 / alpha * log(1 + vV[ ,2] * (exp(-alpha) - 1) / (exp(-alpha * vV[ ,1]) * (1 - vV[ ,2]) + vV[ ,2])))

(mean(vU[, 1] * vU[, 2]) - prod(colMeans(vU))) / sqrt( var(vU[, 1]) * var(vU[, 2]) )

# (c) v

vU.prime <- 1 - vU
(mean(vU.prime[, 1] * vU.prime[, 2]) - prod(colMeans(vU.prime))) / sqrt( var(vU.prime[, 1]) * var(vU.prime[, 2]) )

# (c) vi

X1X2 <- qlnorm(vU, - 0.5 * log(5), sqrt(log(5)))
X1X2.prime <- qlnorm(vU.prime, - 0.5 * log(5), sqrt(log(5)))


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

F1 <- function(x) plnorm(x, - 0.5 * log(5), sqrt(log(5)))
F2 <- function(x) plnorm(x, - 0.5 * log(5), sqrt(log(5)))

rhox1x2 <- function(alpha){
  C.frank <- function(u, v) -1 / alpha * log(1 + (exp(-alpha * u) - 1) * (exp(-alpha * v) - 1) / (exp(-alpha) - 1))
  Fx1x2 <- function(x, y) C.frank(F1(x), F2(y))
  Fx1x2.bar <- function(x, y) 1 - F1(x) - F2(y) + Fx1x2(x, y)
  Ex1x2 <- quad2d(Fx1x2.bar, 0, 100, 0, 100)
  covx1x2 <- Ex1x2 - EX ** 2
  covx1x2 / sqrt(VX ** 2)
}

TrouverAlpha <- function(rho){
  optimize(function(x) abs(rhox1x2(x) - rho), c(0, 30))$minimum
}

(alpha <- TrouverAlpha(0.5))

# (c) iv

nsim <- 1000000
set.seed(2017)
vV <- matrix(runif(nsim * 2), nsim, 2, byrow = TRUE)
vU <- cbind(vV[ ,1], -1 / alpha * log(1 + vV[ ,2] * (exp(-alpha) - 1) / (exp(-alpha * vV[ ,1]) * (1 - vV[ ,2]) + vV[ ,2])))


(mean(vU[, 1] * vU[, 2]) - prod(colMeans(vU))) / sqrt( var(vU[, 1]) * var(vU[, 2]) )

# (c) v

vU.prime <- 1 - vU
(mean(vU.prime[, 1] * vU.prime[, 2]) - prod(colMeans(vU.prime))) / sqrt( var(vU.prime[, 1]) * var(vU.prime[, 2]) )

# (c) vi

X1X2 <- qlnorm(vU, - 0.5 * log(5), sqrt(log(5)))
X1X2.prime <- qlnorm(vU.prime, - 0.5 * log(5), sqrt(log(5)))

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
  C.frank <- function(u, v) -1 / alpha * log(1 + (exp(-alpha * u) - 1) * (exp(-alpha * v) - 1) / (exp(-alpha) - 1))
  C.frank.bar <- function(u, v) u + v - 1 + C.frank(1 - u, 1 - v)
  Fx1x2 <- function(x, y) C.frank.bar(F1(x), F2(y))
  Fx1x2.bar <- function(x, y) 1 - F1(x) - F2(y) + Fx1x2(x, y)
  Ex1x2 <- quad2d(Fx1x2.bar, 0, 100, 0, 100)
  covx1x2 <- Ex1x2 - EX ** 2
  covx1x2 / sqrt(VX ** 2)
}

TrouverAlpha <- function(rho){
  optimize(function(x) abs(rhox1x2(x) - rho), c(0, 30))$minimum
}

(alpha <- TrouverAlpha(0.5))

# (c) iv

nsim <- 1000000
set.seed(2017)
vV <- matrix(runif(nsim * 2), nsim, 2, byrow = TRUE)
vU <- cbind(vV[ ,1], -1 / alpha * log(1 + vV[ ,2] * (exp(-alpha) - 1) / (exp(-alpha * vV[ ,1]) * (1 - vV[ ,2]) + vV[ ,2])))


(mean(vU[, 1] * vU[, 2]) - prod(colMeans(vU))) / sqrt( var(vU[, 1]) * var(vU[, 2]) )

# (c) v

vU.prime <- 1 - vU
(mean(vU.prime[, 1] * vU.prime[, 2]) - prod(colMeans(vU.prime))) / sqrt( var(vU.prime[, 1]) * var(vU.prime[, 2]) )

# (c) vi

X1X2 <- qlnorm(vU, - 0.5 * log(5), sqrt(log(5)))
X1X2.prime <- qlnorm(vU.prime, - 0.5 * log(5), sqrt(log(5)))

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


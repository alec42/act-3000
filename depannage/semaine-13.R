# Code R presente lors du 1er decembre 2017
# Christopher Blier-Wong
# ACT-3000 Theorie du risque

# Plan du depannage -------------------------------------------------------

# Estimation

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

# Maximum de vraisemblance

Densite.EFGM <- function(u1, u2, alpha){
  1 + alpha * (1 - 2 * u1) * (1 - 2 * u2)
}

NegLogVraisemblance <- function(alpha){
  - sum(sapply(1:nrow(vU), function(t) log(Densite.EFGM(vU[t, 1], vU[t, 2], alpha))))
}

alpha.span <- seq(-0.9, 0.9, 0.1)
plot(alpha.span, sapply(alpha.span, NegLogVraisemblance), type = "l")

optimize(NegLogVraisemblance, c(-1, 1))$minimum

# Estimation d'une paire de variables aléatoires

vX <- qexp(vU, 1)

# Methode 1: methode directe

NegLogVraisemblance <- function(theta){
  - sum(sapply(1:nrow(vU), function(t) log(dexp(vX[t, 1], theta[1]) * dexp(vX[t, 2], theta[2]) * Densite.EFGM(pexp(vX[t, 1], theta[1]), pexp(vX[t, 2], theta[2]), theta[3]))))
}

optim(c(2, .5, 0), NegLogVraisemblance)$par

# Methode 2: IFM Joe

lam.ml <- optimize(function(t) - sum(log(dexp(vX[, 1], t))), c(0, 10))$minimum
(lam.ml <- c(lam.ml, optimize(function(t) - sum(log(dexp(vX[, 2], t))), c(0, 10))$minimum))

vU.2 <- sapply(1:2, function(t) pexp(vX[, t], lam.ml[t]))
vU.2 <- pexp(vX, lam.ml) # R est assez intelligent pour faire ca.

NegLogVraisemblance <- function(alpha){
  - sum(sapply(1:nrow(vU.2), function(t) log(Densite.EFGM(vU.2[t, 1], vU.2[t, 2], alpha))))
}

optimize(NegLogVraisemblance, c(-1, 1))$minimum

# Methode 3: Semi-parametrique

vU.3 <- cbind(rank(vX[, 1]) / (nrow(vX) + 1), rank(vX[, 2]) / (nrow(vX) + 1))

NegLogVraisemblance <- function(alpha){
  - sum(sapply(1:nrow(vU.3), function(t) log(Densite.EFGM(vU.3[t, 1], vU.3[t, 2], alpha))))
}

optimize(NegLogVraisemblance, c(-1, 1))$minimum

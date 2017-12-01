# Code R presente lors du 1er decembre 2017
# Christopher Blier-Wong
# ACT-3000 Theorie du risque

# Plan du depannage -------------------------------------------------------

# Suite des copules
# Estimation couples
# Allocation du capital 

# Suite des copules -------------------------------------------------------

# 10.5.3

# 10.5.3

# (d)

alpha <- 1 / 2

Fs <- function(x) alpha * pexp(x / 2, 1) + (1 - alpha) * ifelse((x > log(4)), sqrt(1 - 4 * exp(-x)), 0) + 
  (1 - alpha - 1 + alpha) * pgamma(x, 2)

VaR <- function(kap) {
  return(optimize(function(s) abs(Fs(s) - kap), c(0, 100))$minimum)
}
VaR(0.01)
VaR(0.99)

g <- function(n) {
  o <- vector("numeric", length(n))
  for(i in seq(along = n))
    o[i] <- optimize(VaR, c(0, 80))$minimum
  n * o
}

# pas une question pour R, integration pourrie 
sapply(c(0.01, 0.99), function(t) integrate(g, lower = t, upper = 1)$value / (1 - t)) 

# (e)

integrale_que_je_nai_pas_le_gout_de_faire <- integrate(function(x) x * log(1 - exp(- x)) * exp(- x), 0, 100)$value

alpha <- (1 + integrale_que_je_nai_pas_le_gout_de_faire) / (2 + integrale_que_je_nai_pas_le_gout_de_faire)

Fs <- function(x) alpha * pexp(x / 2, 1) + (1 - alpha) * ifelse((x > (log(4))), sqrt(1 - 4 * exp(-x)), 0) 

par(mfrow=c(1, 2))
curve(Fs, 0, 40)
curve(pgamma(x, 2, 1), 0, 40)
par(mfrow=c(1, 1))

VaR <- function(kap) optimize(function(s) abs(Fs(s) - kap), c(0, 100))$minimum
VaR(0.01)
VaR(0.99)

g <- function(n) {
  o <- vector("numeric", length(n))
  for(i in seq(along = n))
    o[i] <- optimize(VaR, c(0, 100))$minimum
  n * o
}

# pas une question pour R, integration pourrie
sapply(c(0.01, 0.99), function(t) integrate(g, lower = t, upper = 1)$value / (1 - t)) 

# (f)

sig <- 1
mu <- log(10) - sig / 2

integrale_que_je_nai_pas_le_gout_de_faire_2 <- integrate(function(x) x * exp(mu + sig * qnorm(1 - plnorm(x, mu, sig))) * dlnorm(x, mu, sig), 0, 100)$value
alpha <- (exp(2 * mu + sig ** 2) - integrale_que_je_nai_pas_le_gout_de_faire_2) / (exp(2 * mu + 2 * sig ** 2) - integrale_que_je_nai_pas_le_gout_de_faire_2)

Fs <- function(s) alpha * plnorm(s, mu, sig) + (1 - alpha) * 
  ifelse((s > 2 * exp(mu)), pnorm(1 / sig * log( (s + sqrt(s ** 2 - 4 * exp(2 * mu))) / (2 * exp(mu)))) - pnorm(1 / sig * log( (s - sqrt(s ** 2 - 4 * exp(2 * mu))) / (2 * exp(mu)))), 0)

curve(Fs, 0, 100)

VaR <- function(kap) optimize(function(s) abs(Fs(s) - kap), c(0, 100))$minimum
VaR(0.01)
VaR(0.99)

g <- function(n) {
  o <- vector("numeric", length(n))
  for(i in seq(along = n))
    o[i] <- optimize(VaR, c(0, 90))$minimum
  n * o
}

# pas une question pour R, integration pourrie
sapply(c(0.01, 0.99), function(t) integrate(g, lower = t, upper = 1)$value / (1 - t)) 

# 10.5.5

# (d)

nsim <- 1000000
set.seed(2018)

theta <- floor(runif(nsim) * 4) + 1
# theta <- sample(1:4, nsim, replace = TRUE)
alpha <- 5
vV <- matrix(runif(nsim * 3), nsim, 3, byrow = T)
vTheta <- qgamma(vV[, 1], 1 / alpha, 1)
vY <- sapply(1:2, function(t) qexp(vV[, t + 1], vTheta))
U <- (1 + vY) ** (-1 / alpha)

W <- matrix(numeric(), nsim, 2)
for (i in 1:nsim) {
   if (theta[i] == 1) {
     W[i, ] <- c(U[i, 1], U[i, 2])
   } else if (theta[i] == 2) {
     W[i, ] <- c(U[i, 1], 1 - U[i, 2])
   } else if  (theta[i] == 3) {
     W[i, ] <- c(U[i, 1], 1 - U[i, 2])
   } else {
     W[i, ] <- c(1 - U[i, 2], 1 - U[i, 2])
   }
}

# (i)

vT <- rowSums(W)

plot.ecdf(vT[1:100000])

toutes_les_questions_sauf_i <- function(quantile_function) {
  X <- quantile_function(W)
  S <- rowSums(X)
  plot.ecdf(S[1:10000])
  S.ind <- quantile_function(runif(nsim)) + quantile_function(runif(nsim))
  plot.ecdf(S.ind[1:10000])
  print(paste0("VaR_", kappas, "(S) = ", sort(S)[nsim * kappas]))
  print(paste0("TVaR_", kappas, "(S) = ", sapply(kappas, function(t) mean(sort(S)[(t * nsim):nsim]))))
}

q1 <- function(x) qnorm(x, 0, 1)
q2 <- function(x) qlnorm(x, -1/2, 1)
q3 <- function(x) qexp(x, 1)
q4 <- function(x) 1.5 * ((1 - x) ** (- 1 / 2.5) - 1)

kappas <- c(1/1000, 1/100, 99/100, 999/1000)

par(mfrow = c(1, 2))
toutes_les_questions_sauf_i(q1)
toutes_les_questions_sauf_i(q2)
toutes_les_questions_sauf_i(q3)
toutes_les_questions_sauf_i(q4)
par(mfrow = c(1, 1))

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

# Allocation du capital  --------------------------------------------------

# 11.4.1

# b)

aa <- 2 ** 5

fb <- c(0, 1, rep(0, 30))
phi1 <- fft(fb)
phin <- (0.9 + 0.1 * phi1) * phi1 * (0.06 + 0.02 * (0.8 + 0.2 * phi1) ** 2)
E1 <- Re(fft(phin, inverse = TRUE)) / aa
round(E1[1:5], 5)

fb <- c(0, 1, rep(0, 30))
phi1 <- fft(fb)
phin <- (0.8 + 0.2 * phi1) * phi1 * (0.08 + 0.04 * (0.9 + 0.1 * phi1) ** 2)
E2 <- Re(fft(phin, inverse = TRUE)) / aa
round(E2[1:5], 5)

(round(E1[1:5], 5) + round(E2[1:5], 5)) / (0:4)


# 11.4.2

fb <- c(0.8, 0.16, 0.04)
fb <- c(fb, rep(0, aa - length(fb)))

phib <- fft(fb)
phin <- phib ** 12
fn <- Re(fft(phin, inverse = TRUE)) / aa
sum(fb)
fn[6]

fb <- c(0, 1, rep(0, 30))
phi1 <- fft(fb)
phin <- (0.8 + 0.16 * phi1 + 0.04 * phi1 ** 2) ** 12
fn <- Re(fft(phin, inverse = TRUE)) / aa
sum(fb)
fn[6]

# b)

fb <- c(0, 1, rep(0, 30))
phi1 <- fft(fb)
phin <- phi1 * 12 * (0.1 + 0.04 * phi1) * (0.8 + 0.16 * phi1 + 0.04 * phi1 ** 2) ** 11
E1 <- Re(fft(phin, inverse = TRUE)) / aa
E1[6]

fb <- c(0, 1, rep(0, 30))
phi1 <- fft(fb)
phin <- phi1 * 12 * (0.06 + 0.04 * phi1) * (0.8 + 0.16 * phi1 + 0.04 * phi1 ** 2) ** 11
E2 <- Re(fft(phin, inverse = TRUE)) / aa
E2[6]

E1[2] + E2[2]
fn[2]

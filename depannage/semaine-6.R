
# Document préparatoire pour l'examen informatique ------------------------

# Christopher Blier-Wong
# 2017-10-20

# Dans ce document, je vais faire quelques exemples de codes. Il ne faud 
# pas apprendre par coeur, mais comprendre ce qu'il se passe. 

# Rappel de l'aide en R ---------------------------------------------------

?optimize
help(optimize)

# click + f1 pour l'aide
# click + f2 pour la fonction

# Autres rappels ----------------------------------------------------------

c(1, 2, 3) # creer un vecteur
1:3 # cree une suite de 1 a 3
# & est le "et" logique
# | est le "ou" logique
# % est modulo
# %/% est la division entiere
# F = FALSE = 0
# T = TRUE = 1

Y1 <- c(0, 0.1, 0.1, 0.4, 0.5, 0.6, 0.7, 0.7, 0.7, 0.9, 1, 1.3, 4.5, 5.6)

(Y1 > 1) # retourne vecteur de vrai/faux
sum(Y1 > 1) # retourne le nombre d'elements plus grand que 1
mean(Y1 > 1) # retourne la probabilite empirique d'etre plus grand que 1
Y1[Y1 > 1]
mean((Y1 > 0.5) & (Y1 < 1)) # retourne la probabilite empirique d'etre entre 0.5 et 1
mean((Y1 > 0.5) * (Y1 < 1)) # meme chose

# Rappel simulation -------------------------------------------------------

# Souvent en examen, on demande de simuler plusieurs variables aleatoires en
# meme temps. Voici un rappel. 

# Methode lente

set.seed(2017)
m <- 1000000
u <- runif(m*2)

U <- matrix(numeric(), m, 2)

for(i in 1:m){
  U[i, 1] <- u[2*(i - 1) + 1]
  U[i, 2] <- u[2*(i - 1) + 2]
}

U[1234, ]

# Methode rapide

set.seed(2017)
U <- matrix(runif(2*m), m, 2, byrow = TRUE)
U[1234, ]


# Rappel GNPA -------------------------------------------------------------

# GNPA de l'examen partiel informatique A16

a <- 41358
m <- 2**31 - 1
x0 <- 20150418
mm <- 10000


x <- c(x0, numeric(3*mm))

for(i in 1 + 1:(3*mm)){
  x[i] <- (a*x[i-1]) %% m
}

u <- x[-1]/m
U <- matrix(u, mm, 3, byrow = TRUE)
U[c(1, 1000, 10000), ]



# Simulation --------------------------------------------------------------

# Si les v.a. suivent la meme loi, on peut tout simplement inverser

X <- qexp(U, 0.2)
X[1, ]


# Sinon, on doit traiter chaque colonnes individuellement
# Examen informatique partiel A16, question 1

Y1 <- qgamma(U[, 1], 2, 1)
Y2 <- qgamma(U[, 2], 1/2, 1)
Y3 <- qgamma(U[, 3], 1, 1)

X <- cbind(100*(Y1 + Y3), 200*(Y2 + Y3))
X[2, ] # verif ok
X[c(3, 4)] # reponse

# Approximations avec simulation ------------------------------------------

# Avec des simulations, on peut tout approximer. 
# suite Examen informatique partiel A16, question 1

#Fx1x2(500, 500)

mean((X[, 1] < 500) & (X[, 2] < 500))

# pmax et pmin prennent le minimum par colonne. 
# Utile pour trouver des esperances de max ou min
# Pourrait aussi etre utile pour trouver les bornes
# de Frechet

mean(pmax(X[, 1] - 500, 0) * (X[, 2] > 500))

S <- rowSums(X)

# VaR : on prend le kappa * m ieme element du vecteur ordonne
s.sorted <- sort(S)
(VaRS <- s.sorted[0.99*mm])

# TVaR : on prend la moyenne des valeurs de S plus grand que VaR
mean(S[S > VaRS])
mean(s.sorted[(0.99*mm+1):mm])



# Rappel optimisation numérique -------------------------------------------

# Pour trouver le minimum d'une fonction, on prend la fonction optimize.

help(optimize)

curve((x-3)^2, 0, 5)

(res_optim <- optimize(function(x) (x-3)^2, c(0, 5)))
res_optim$minimum
res_optim$objective


# Inversion numerique -----------------------------------------------------

# On cherche la VaR pour une distribution Gamma avec kappa = 0.9

curve(pgamma(x, 1, 1), 0, 10, xaxs = "i", yaxs = "i")
curve(pgamma(x, 1, 1) - 0.9, 0, 10, xaxs = "i", yaxs = "i")
abline(h = 0)
curve(abs(pgamma(x, 1, 1) - 0.9), 0, 10, xaxs = "i", yaxs = "i")
axis(1, qgamma(0.9, 1, 1), expression(VaR[kappa](X)), col = 1)

axis(1, qgamma(0.9, 1, 1), expression(VaR[kappa](X)), col = 1)

F_to_minimize <- function(x) abs(pgamma(x, 1, 1) - 0.9)

(VaRX <- optimize(F_to_minimize, c(0, 5))$minimum)
qgamma(0.9, 1, 1) # verification 1
pgamma(VaRX, 1, 1) # verification 2


# On cherche la fontion de repartition pour des v.a. 
# antimonotones exponentielles 

VaRtemp <- function(u) qexp(u, 0.2) + qexp(1 - u, 0.2)
curve(VaRtemp, 0, 1)

VaRS <- function(u) qexp((1 - u)/2, 0.2) + qexp((1 + u)/2, 0.2)
curve(VaRS, 0, 1, xaxs = "i")
axis(2, VaRS(0), expression(VaR[0](S)), col = 1)

F_to_minimize <- function(x) abs(VaRS(x) - 15)
optimize(F_to_minimize, c(0, 1))$minimum

# Attention aux erreurs de minimisation
optimize(F_to_minimize, c(0, 1))$minimum
optimize(F_to_minimize, c(0.8, 1))$minimum
optimize(F_to_minimize, c(0.89, 0.9))$minimum
(Fx <- optimize(F_to_minimize, c(0.89, 0.9), tol = .Machine$double.eps)$minimum)
VaRS(Fx) # verification

# Retour sur la simulation ------------------------------------------------

# Pour des lois multivariees, on peut inverser numeriquement au lieu de
# inverser analytiquement.
# Exemple EFGM beta1 = 1, beta2 = 2, theta = 0.5

m <- 10000

U <- matrix(runif(2*m), m, 2, byrow = TRUE)


Fx1x2 <- function(x1, x2) (1 - exp(-x1))*(1 - exp(-2*x2)) +
  0.5*(1 - exp(-x1))*(1 - exp(-2*x2))*exp(-x1)*exp(-2*x2)

X <- cbind(qexp(U[, 1], 1), numeric(mm))

for(i in 1:mm){
  X[i, 2] <- optimize(function(t) abs(Fx1x2(X[i, 1], t) - U[i, 2]), c(0, 10))$minimum
}

# Domaines d'optimisation -------------------------------------------------

# Il faut faire attention au domaine d'optimisation. 

# Pour une VaR, on peut optimiser kappa sur (0, 1)
# Pour une fonction de repartition, on peut optimiser x sur le support de x
# Pour une fgm, on peut optimiser sur le domaine de t (voir annexe)

# Quand on optimise, on peut avoir une reponse plus precise en
# 1. Raptisser l'intervalle d'optimisation
# 2. Reduire le parametre de tolerance


# Algorithmes recrusifs ---------------------------------------------------

## DePril

recur.nrisks<-function(ff,nn=5,smax=100)
{
  # convolution de n fns de masses de probabilité avec
  # elle-meme
  # premier algorihtme de DePril
  ll<-length(ff)
  ffs<-ff[1]^nn
  ff<-c(ff,rep(0,smax-ll+1))
  27
  for (i in 1 :smax)
  {
    j<-i+1
    ffs<-c(ffs,(1/ff[1])*sum(ff[2 :j]*ffs[i :1]*((nn+1)*(1 :i)/i-1)))
  }
  return(ffs)
}

fb <- dpois(0:100, 2)

# smax >= longueur du vecteur fb
recur.nrisks(fb, nn=5, smax=5)
recur.nrisks(fb, nn=5, smax=99)
recur.nrisks(fb, nn=5, smax=100)[1:5]

dpois(0:4, 2 * 5)

# Panjer

panjer.poisson<-function(lam,ff,smax)
{
  aa<-0
  bb<-lam
  ll<-length(ff)
  ffs<-exp(lam*(ff[1]-1))
  ff<-c(ff,rep(0,smax-ll+1))
  for (i in 1 :smax)
  {
    j<-i+1
    ffs<-c(ffs,(1/(1-aa*ff[1]))*sum(ff[2 :j]*ffs[i :1]*(bb*(1 :i)/i+aa)))
  }
  return(ffs)
}

fb <- c(0, 1)
panjer.poisson(10, fb, smax=4)
dpois(0:4, 2 * 5)

fb <- dbinom(0:10, 10, 0.2)
panjer.poisson(10, fb, smax=10)

# FFT

# Vous pouvez prendre les codes en annexe, mais le code est simple

# Calculer les valeurs pour k = 0, 1, ..., 16383
# Il vous donne la valeur de m prendre. On peut prendre tous entiers
# plus grands que m = 14

(m <- log(16384)/log(2))
a <- 2^m # ou a <- 16384 directement, si vous savez que 16384 est 2^x

fb <- dbinom(0:10, 10, 0.2)
fb <- c(fb, numeric(a - length(fb)))
length(fb) # verif ok
phib <- fft(fb)


# Somme de deux v.a.
phix <- phib^2
fx <- Re(fft(phix, inverse=TRUE))/a
fx[1:5]
dbinom(0:4, 20, 0.2) # verif ok

# Somme aleatoire

phix <- exp(4*(phib - 1)) # Pm(phib(t))
fx <- Re(fft(phix, inverse=TRUE))/a
fx[1:5]
panjer.poisson(4, fb[1:11], smax=10)[1:5] # Panjer = fft

# Loi multivariee

Pm1m2 <- function(t1, t2) (0.8 + 0.1*t1 + 0.02*t2 + 0.08*t1*t2)**12

fb1 <- dbinom(0:10, 10, 0.2)
fb2 <- c(0, dgeom(0:100, 0.2))

fb1 <- c(fb1, numeric(a - length(fb1)))
fb2 <- c(fb2, numeric(a - length(fb2)))

phib1 <- fft(fb1, inverse=FALSE)
phib2 <- fft(fb2, inverse=FALSE)

phix <- Pm1m2(phib1, phib2) 
fx <- Re(fft(phix, inverse=TRUE))/a
fx[1:5]

# Verifications et rappels ------------------------------------------------

phix <- exp(4*(phib - 1)) # Pm(phib(t))
fx <- Re(fft(phix, inverse=TRUE))/a

# Vu que fx(0) est dans le premier element du vecteur, on doit faire + 1 
# pour trouver fx(k)

fx[1] # fx(0)
fx[2] # fx(1)
fx[1 + c(5, 10, 15, 20)] # fx(5), fx(10), fx(15), fx(20)

# DePril, Panjer, fft. Assurez-vous que la somme des elements donne 1
sum(fx)

# Si la moyenne est disponible, vous pouvez l'utiliser pour valider.

sum(fx*(0:(a - 1)))
4*10*0.2 # Poisson(4) et Binom(10, 0.2)

# La fonction de repartition est la somme cumulative des fx

Fx <- cumsum(fx)

Fx[1] # Fx(0)
Fx[2] # Fx(1)
Fx[1 + c(5, 10, 15, 20)] # Fx(5), Fx(10), Fx(15), Fx(20)

# La VaR est la plus petite valeur de X qui satisfait Fx(x) > kappa
# On retire 1 parce que Fx(0) = Fx[1]

min(which(Fx > 0.99)) - 1

sapply(c(0.9, 0.99, 0.999), function(t) min(which(Fx > t)) - 1)

# On peut tout calculer car on a la fonction de densite

sl <- pmax(0:(a-1) - 10, 0)
sl[1:20]

# E[max(X-10, 0)]

sum(sl*fx)

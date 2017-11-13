# Ce code est le solutionnaire pour l'examen partiel informatique A17
# Code developpe par les auxiliaires Christopher Blier-Wong et 
# Ihsan Chaoubi. 

# Question 1 --------------------------------------------------------------

FX1 <- function(x) exp(-(x/1000)^(-2))
FX2 <- function(x) 1 + (log(1 - 0.5*exp(-x/2000))/log(2))

VaRX1 <- function(kap) sqrt(-(1000^2)/log(kap))
VaRX2 <- function(kap) -2000 * log(2 - 2^kap)
VaRS <- function(kap) VaRX1(kap) + VaRX2(kap)

FX1(VaRX1(0.9)) # check ok
FX2(VaRX2(0.9)) # check ok

VaRcx <- function(VaR, kap) 0.25*VaR(1 - 0.25*(1-kap)) + 0.5*VaR(1 - 0.5*(1-kap)) + 0.25*VaR(1 - 0.75*(1-kap))

VaRcx(VaRX1, 0.9) # reponse
VaRcx(VaRX2, 0.9) # reponse

VaRcx(VaRX1, 0.9) + VaRcx(VaRX2, 0.9) # reponse
VaRcx(VaRS, 0.9) # reponse alternative

(Fs5000 <- optimize(function(t) abs(VaRX1(t) + VaRX2(t) - 5000), c(0.81, 0.812))$minimum)
VaRS(Fs5000) # check ok

# Question 2 --------------------------------------------------------------

params <- c(log(100) - 0.5, log(200) - 0.32, 1, 0.8, 0.5)

m <- 1000000

set.seed(2017)
U <- matrix(runif(m*2), m, 2, byrow = TRUE)
Z <- qnorm(U)

Y1 <- params[1] + params[3]*Z[, 1]
Y2 <- params[2] + params[4]*(params[5]*Z[, 1] + sqrt(1 - params[5]^2)*Z[, 2])

X <- exp(cbind(Y1, Y2))
S <- rowSums(X)

U[c(1, m), ] # verif
Z[c(1, m), ] # verif

U[2, ] # reponse i
Z[2, ] # reponse i
X[2, ] # reponse i

X1.sort <- sort(X[, 1])
X2.sort <- sort(X[, 2])

# reponses ii
X1.sort[m*0.99] # VaR X1
X2.sort[m*0.99] # VaR X2
0.25*X1.sort[m*(1-1/4*(1-.99))] + 0.5*X1.sort[m*(1-2/4*(1-.99))] + 0.25*X1.sort[m*(1-3/4*(1-.99))] # VaR convexe X1
0.25*X2.sort[m*(1-1/4*(1-.99))] + 0.5*X2.sort[m*(1-2/4*(1-.99))] + 0.25*X2.sort[m*(1-3/4*(1-.99))] # VaR convexe X2

mean(X1.sort[(0.99*m + 1):m]) # TVaR X1
mean(X2.sort[(0.99*m + 1):m]) # TVaR X2

# reponses iii

S.sort <- sort(S)

S.sort[m*0.99] # VaR S
0.25*S.sort[m*(1-1/4*(1-.99))] + 0.5*S.sort[m*(1-2/4*(1-.99))] + 0.25*S.sort[m*(1-3/4*(1-.99))] # VaR convexe S
mean(S.sort[(0.99*m + 1):m]) # TVaR S

# reponse iv

0.25*X1.sort[m*(1-1/4)] + 0.5*X1.sort[m*(1-2/4)] + 0.25*X1.sort[m*(1-3/4)] # VaR convexe X1
0.25*X2.sort[m*(1-1/4)] + 0.5*X2.sort[m*(1-2/4)] + 0.25*X2.sort[m*(1-3/4)] # VaR convexe X2
0.25*S.sort[m*(1-1/4)] + 0.5*S.sort[m*(1-2/4)] + 0.25*S.sort[m*(1-3/4)] # VaR convexe S

0.25*X1.sort[m*(1-1/4)] + 0.5*X1.sort[m*(1-2/4)] + 0.25*X1.sort[m*(1-3/4)] + 
  0.25*X2.sort[m*(1-1/4)] + 0.5*X2.sort[m*(1-2/4)] + 0.25*X2.sort[m*(1-3/4)] 
# VaRcxS > VaRcxX1 + VaRcxX2 alors pas sous-additive

# reponse v

mean(X1.sort[(0.999*m + 1):m]) + mean(X2.sort[(0.999*m + 1):m]) - mean(S.sort[(0.999*m + 1):m])

# Question 3 --------------------------------------------------------------

params <- c(10, 0.2, 10, 0.3, 0.8, 0.3, 0.2)

fb1 <- dbinom(0:params[1], params[1], params[2])
fb2 <- dbinom(0:params[3], params[3], params[4])
Pm12 <- function(t1, t2) exp(params[5]*(t1 - 1))*exp(params[6]*(t2 - 1))*exp(params[7]*(t1*t2 - 1))

m <- log(4096)/log(2)
aa <- 2**m

phib1 <- fft(c(fb1, numeric(aa - length(fb1))))
phib2 <- fft(c(fb2, numeric(aa - length(fb2))))

phis <- Pm12(phib1, phib2)

fs <- Re(fft(phis, inverse = TRUE)/aa)

fs[1]
sum(fs*0:4095) # moyenne ok

fs[1 + c(2, 3, 4)] # verif ok
fs[1 + c(0, 1, 5, 10)] # reponse

Fs <- cumsum(fs)

min(which(Fs > 0.9999)) - 1 # verif ok
sapply(c(0.9, 0.99, 0.999), function(t) min(which(Fs > t)) - 1) # reponse

# Question 4 --------------------------------------------------------------

params <- c(1, 0.5, 1)


M1 <- function(t, bet) (bet/(bet - t))
Ms <- function(t){
  (1 + params[3])*M1(t, params[1])*M1(t, params[2]) - 
    params[3]*M1(t, 2*params[1])*M1(t, params[2]) -
    params[3]*M1(t, params[1])*M1(t, 2*params[2]) + 
    params[3]*M1(t, 2*params[1])*M1(t, 2*params[2])
} 

phi <- function(t, kap) 1/t*log(Ms(t)/(1 - kap))

xx <- seq(0.1, 0.49, length.out = 1000)
plot(xx, sapply(xx, function(t) phi(t, 0.99)), type = "l") # fonction convexe

min_res <- optimize(function(t) phi(t, 0.99), lower = 0, upper = .5, tol = .Machine$double.eps)

min_res$minimum # parametre tk
min_res$objective # eVaR

# Question 5 --------------------------------------------------------------

Fx <- function(x) 1 - (2/(2 + x))^1.2
VaRX <- function(kap) 2*((1 - kap)^(-1/1.2) - 1)
VaRS <- function(kap) VaRX(0.5 - kap/2) + VaRX(0.5 + kap/2)

optimize(function(t) abs(VaRS(t) - 300), c(0.99, 1))$minimum
optimize(function(t) abs(VaRS(t) - 100), c(0.98, 0.99))$minimum

# Question 6 --------------------------------------------------------------

# a) 

(scale1 <- -1/3*log(0.01))
(scale2 <- -1/2*log(0.01))

# b) 

1/scale1 # moyenne W1
1/scale2 # moyenne W2

# d) 

(scalen <- scale1 + scale2)

fc <- c(0, 0.4*dgeom(0:100, 0.2) + 0.6*dgeom(0:100, 0.4))

# e) 

1/scale2 # moyenne Wn

# f)

m <- 14
aa <- 2^m

fc <- c(fc, numeric(aa - length(fc)))
phic <- fft(fc, inverse = FALSE)
phin <- exp(4*scalen*(phic - 1))
fs <- Re(fft(phin, inverse = TRUE)/aa)

fs[1 + c(50, 60, 70)] # reponse densite


# ici on prend l'indice k + 1 + 1 parce que fs[1] est f(0)
# et on prend un de plus parce que c'est plus grand et non egal
sum(fs[(100 + 1 + 1):aa]) # verification ok
sum(fs[(105 + 1 + 1):aa]) # reponse

# g) 
Fs <- cumsum(fs)
# moins un parce que le premier element est pour f(0)
min(which(Fs > 0.99)) - 1 # verification ok
min(which(Fs > 0.9)) - 1 # reponse

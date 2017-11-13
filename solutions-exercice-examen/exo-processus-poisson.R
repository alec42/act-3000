# Reponses aux questions sur les processus de poisson

# Question 3 --------------------------------------------------------------

(scale <- -1/3*log(0.01))
1/scale # Esperance inter-sinistre
dpois(2, 2*scale) # P(N(2) = 2)

fb <- c(0, dgeom(0:100, 0.5))

fs <- panjer.poisson(4*scale, fb, 101)

Fs <- cumsum(fs)
1 - Fs[1 + 6] # reponse desiree
1 - Fs[1 + 20] # reponse dans solutionnaire


# Question 4 --------------------------------------------------------------

scale = 1

Y1 <- c(1.587, 0.307, 0.611, 1.474, 1.422)
Y2 <- c(0.589, 0.777, 4.185, 4.922, 5.424)

# T2 ~ Gamma(2, 1)
1 - pgamma(3, 2, 1) # P(T2 > 3)


# T3 ~ Gamma(3, 1)
3/1 # Esperance
qgamma(0.5, 3, 1) # Mediane


t <- c(1:10)/2
T1 <- cumsum(Y1)
sapply(t, function(j) sum(T1 < j))

T2 <- cumsum(Y2)
sapply(t, function(j) sum(T2 < j))

# Question 5 --------------------------------------------------------------

a <- 16384
scale <- 4.5

fb <- c(0, 4, 0.5)/scale
fb <- c(fb, numeric(a - length(fb)))

phib <- fft(fb, inverse=FALSE)
phis <- exp(scale*1*(phib - 1))
fx <- Re(fft(phis, inverse = TRUE))/a
fx[1:11]

phis <- exp(scale*1.75*(phib - 1))
fx <- Re(fft(phis, inverse = TRUE))/a
fx[1:11]

# Question 6 --------------------------------------------------------------

a <- 16384
scalen <- 0.15
fb <- c(0, 0.05, 0.04, 0.03, 0.02, 0.01)/0.15

fs <- panjer.poisson(1*scalen, fb, 7)
fs[1:4]

# Construction d'un code fft maison
# Base sur un code dans les annexes pour ACT-3000 A16

fftmaison <- function(ff){
  nn <- length(ff)
  vk <- 1:nn - 1
  ft <- sapply(vk, function(t) sum(ff*exp(2i*pi*t/nn*vk)))
  return(ft)
}

fftinvmaison <- function(ft){
  nn <- length(ft)
  vk <- 1:nn - 1
  ff <- sapply(vk, function(t) sum(ft*exp(-2i*pi*t/nn*vk)))
  return(Re(ff)/nn)
}

# Exemple avec une loi binomiale
n1 = 3; p = 0.1; n = 10

fbin <- dbinom(0:n1, n1, p)
fbin <- c(fbin, numeric(2**7 - length(fbin)))

phibin <- fftmaison(fbin)
phibin2 <- phibin^n
fbin2 <- fftinvmaison(phibin2)

fbin2[1:10]
dbinom(0:9, n1*n, p)

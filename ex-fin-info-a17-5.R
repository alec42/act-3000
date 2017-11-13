# ACT-3000 Examen info final numero 5
# Christopher Blier-Wong
# 12 octobre 2017

# a) 
scale1 <- -1/3*log(0.01)
scale2 <- -1/2*log(0.01)

# b) 

1/scale1 # moyenne W1
1/scale2 # moyenne W2

# d) 

scalen <- scale1 + scale2

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


# ici on prend l'indice k + 1 + 1 parce que fs[1] est f(0)
# et on prend un de plus parce que c'est plus grand et non egal
sum(fs[(100 + 1 + 1):aa]) # verification ok
sum(fs[(105 + 1 + 1):aa]) # reponse


# g) 
Fs <- cumsum(fs)
# moins un parce que le premier element est pour f(0)
min(which(Fs > 0.99)) - 1 # verification ok
min(which(Fs > 0.95)) - 1 # reponse

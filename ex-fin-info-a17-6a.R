# ACT-3000 Examen info final numero 5
# Christopher Blier-Wong
# 12 octobre 2017

m <- 5
n <- 12
aa <- 2**m

fk <- c(0.8, 0.16, 0.04)
fk <- c(fk, numeric(aa - length(fk)))

phik <- fft(fk, inverse=FALSE)
phin <- phik^n
fn <- Re(fft(phin, invers=TRUE)/aa)

# ici on fait k + 1 car f[1] est f(0)
fn[6] # verification ok
fn[c(3, 8) + 1]

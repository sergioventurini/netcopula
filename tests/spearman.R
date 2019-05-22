n <- 50
K <- 20000

sp <- mu1 <- mu2 <- 0
for (k in 1:K) {
	x <- rnorm(n)
	u1 <- pnorm(x)
	mu1 <- mu1 + sum(u1)
	x <- rnorm(n)
	u2 <- pnorm(x)
	mu2 <- mu2 + sum(u2)
	sp <- sp + sum(u1*u2)
}
12*sp/(K*n) - 3
mu1/(K*n)
mu2/(K*n)

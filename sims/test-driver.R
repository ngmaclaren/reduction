## Working with the adjustments to the mutualistic species model.
## With Ds in the range 0.001 to 1.0, there is no bifurcation with C = 1 and D.tilde = N/10
## I recover the bifurcation with D.tilde = N/20, but not all of it. So, I need to mess with that parameter across the different networks.

library(parallel)
library(igraph)
library(deSolve)
library(latex2exp, lib.loc = "/user/neilmacl/rlocal/")
source("../src/functions.R")
ncores <- detectCores() - 1

network <- "email"
dynamics <- "mutualistic"

load(paste0("../data/", network, ".rda"))
load(paste0("../data/fullstate-", network, ".rda"))

g <- get(network)
N <- vcount(g)
A <- as_adj(g, "both", sparse = FALSE)
k <- degree(g)
## mutualistic_parms$Ds <- seq(0, 3, length.out = lout)
## mutualistic_parms$D.tilde <- N/10

## Y.mutualistic <- with(
##     mutualistic_parms,
##     solve_mutualistic(
##         x = rep(x.init, N), B = B, K = K, C = C, Ds = Ds, D.tilde = D.tilde, E = E, H = H, A = A
##     )
## )
## Y <- Y.mutualistic
Y <- fullstate[[dynamics]]
y <- rowMeans(Y)
n <- floor(log(N))

soln <- with(mutualistic_parms, {
    experiment(n, y, Y, Ds, optimize_weights = FALSE, trace = FALSE)
})

pdf("../img/test.pdf")
with(mutualistic_parms, {
    matplot(Ds, Y, type = "l", lty = 1, lwd = .5, col = adjustcolor(1, alpha.f = .25),
            xlab = TeX(r"(D)", italic = TRUE), ylab = TeX(r"($x^*$)", italic = TRUE),
            cex.lab = 1.75, cex.axis = 1.75)
    lines(Ds, y, lty = 1, lwd = 6, col = 1)
    matlines(Ds, Y[, soln$vs], lty = 1, lwd = 4, col = adjustcolor(2, alpha.f = .5))
    lines(Ds, rowMeans(Y[, soln$vs]), lty = 1, lwd = 6, col = 2)
})
dev.off()

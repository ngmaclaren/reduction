library(parallel)
library(igraph)
library(deSolve)
source("../src/functions.R")
ncores <- detectCores() - 1 # backed off from max to allow for overflow

args <- commandArgs(trailingOnly = TRUE)
network <- args[1]
print(network)
fullstatefile <- paste0("../data/fullstate-", network, ".rda")

load(paste0("../data/", network, ".rda"))
g <- get(network)
N <- vcount(g)
A <- as_adj(g, "both", sparse = FALSE)
k <- degree(g)
mutualistic_parms$D <- N/10

print("solving doublewell...")
Y.doublewell <- with(
    doublewell_parms,
    solve_doublewell(x = rep(x.init, N), r = r, Ds = Ds, A = A)
)

print("solving SIS...")
Y.SIS <- with(
    SIS_parms,
    solve_SIS(x = rep(x.init, N), μ = μ, Ds = Ds, A = A)
)

print("solving gene regulatory...")
Y.genereg <- with(
    genereg_parms,
    solve_genereg(x = rep(x.init, N), B = B, f = f, h = h, Ds = Ds, A = A)
)

print("solving mutualistic species...")
Y.mutualistic <- with(
    mutualistic_parms,
    solve_mutualistic(x = rep(x.init, N), B = B, K = K, Cs = Cs, D = D, E = E, H = H, A = A)
)

## print("solving Wilson-Cowan...")
## Y.wilsoncowan <- with(
##     wilsoncowan_parms,
##     solve_wilsoncowan(E = rep(E.init, N), I = rep(I.init, N), c5s = c5s, A = A, N = N)
## )

fullstate <- list(
    dw = Y.doublewell,
    SIS = Y.SIS,
    genereg = Y.genereg,
    mutualistic = Y.mutualistic#,
    ##wilsoncowan = Y.wilsoncowan
)

attr(fullstate, "graph") <- g
save(fullstate, file = fullstatefile)

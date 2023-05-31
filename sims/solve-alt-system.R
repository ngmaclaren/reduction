library(parallel)
library(igraph)
library(deSolve)
source("../src/functions.R")
ncores <- detectCores() - 1 # backed off from max to allow for overflow

args <- commandArgs(trailingOnly = TRUE)
network <- args[1]
print(network)
fullstatefile <- paste0("../data/fullstate-alt-", network, ".rda")

load(paste0("../data/", network, ".rda"))
g <- get(network)
N <- vcount(g)
A <- as_adj(g, "both", sparse = FALSE)
k <- degree(g)

print("solving doublewell...")
dw1 <- with(
    doublewell_parms,
    solve_doublewell(x = rep(x.init, N), r = c(1, 2, 5), Ds = Ds, A = A)
)
dw2 <- with(
    doublewell_parms,
    solve_doublewell(x = rep(x.init, N), r = c(1, 4, 5), Ds = Ds, A = A)
)

print("solving SIS...")
sis1 <- with(
    SIS_parms,
    solve_SIS(x = rep(x.init, N), μ = .5*μ, Ds = Ds, A = A)
)
sis2 <- with(
    SIS_parms,
    solve_SIS(x = rep(x.init, N), μ = 2*μ, Ds = Ds, A = A)
)

print("solving gene regulatory...")
gr1 <- with(
    genereg_parms,
    solve_genereg(x = rep(x.init, N), B = .5*B, f = f, h = h, Ds = Ds, A = A)
)
gr2 <- with(
    genereg_parms,
    solve_genereg(x = rep(x.init, N), B = 2*B, f = f, h = h, Ds = Ds, A = A)
)

print("solving mutualistic species...")
ms1 <- with(
    mutualistic_parms,
    solve_mutualistic(x = rep(x.init, N), B = B, K = .5*K, Cs = Cs, D = D, E = E, H = H, A = A)
)
ms2 <- with(
    mutualistic_parms,
    solve_mutualistic(x = rep(x.init, N), B = B, K = 2*K, Cs = Cs, D = D, E = E, H = H, A = A)
)

print("solving Wilson-Cowan...")
wc1 <- with(
    wilsoncowan_parms, {
        wilsoncowan_parms$c1 <- wilsoncowan_parms$c1/2
        wilsoncowan_parms$c3 <- wilsoncowan_parms$c3/2
        solve_wilsoncowan(E = rep(E.init, N), I = rep(I.init, N), c5s = c5s, A = A, N = N)
    }
)
wc2 <- with(
    wilsoncowan_parms, {
        wilsoncowan_parms$c2 <- wilsoncowan_parms$c2/2
        wilsoncowan_parms$c4 <- wilsoncowan_parms$c4/2
        solve_wilsoncowan(E = rep(E.init, N), I = rep(I.init, N), c5s = c5s, A = A, N = N)
    }
)

fullstate_alt <- list(
    dw1 = dw1, dw2 = dw2,
    sis1 = sis1, sis2 = sis2,
    gr1 = gr1, gr2 = gr2,
    ms1 = ms1, ms2 = ms2,
    wc1 = wc1, wc2 = wc2
)

attr(fullstate_alt, "graph") <- g
save(fullstate_alt, file = fullstatefile)

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
## mutualistic_parms$D <- N/10

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
    solve_genereg(x = rep(x.init, N), B = B, f = 2*f, h = h, Ds = Ds, A = A) # .5*f involves neg
)
gr2 <- with(
    genereg_parms,
    solve_genereg(x = rep(x.init, N), B = B, f = f, h = 2*h, Ds = Ds, A = A)
)

print("solving mutualistic species...")
ms1 <- with(
    mutualistic_parms,
    solve_mutualistic(
        x = rep(x.init, N), B = B, K = K - 2, C = C, Ds = Ds, D.tilde = D.tilde, E = E, H = H, A = A
    )
)
ms2 <- with(
    mutualistic_parms,
    solve_mutualistic(
        x = rep(x.init, N), B = B, K = K + 2, C = C, Ds = Ds, D.tilde = D.tilde, E = E, H = H, A = A
    )
)

fullstate_alt <- list(
    dw1 = dw1, dw2 = dw2,
    sis1 = sis1, sis2 = sis2,
    gr1 = gr1, gr2 = gr2,
    ms1 = ms1, ms2 = ms2
)

attr(fullstate_alt, "graph") <- g
save(fullstate_alt, file = fullstatefile)

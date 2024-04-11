library(optNS)
options(digits = 3)

load("dolphin-demo.RData")

                                        # ε for dolphin network n \in 1, 2, 3, 4
round(get_error(solns), 3)

                                        # Confirm that we find the exact solution with n = 1 or n = 2
find_exact <- function(n, ...) {
    VS <- combn(N, n)
    errors <- apply(VS, 2, obj_fn, y, Y)

    idx <- which.min(errors)
    vs <- VS[, idx]
    error <- errors[idx]
    ks <- k[vs]

    ## if(optimize_weights) ws <- quadoptm(vs, y, Y)

    ## if(optimize_weights) {
    ##     return(list(vs = vs, error = error, ks = ks, ws = ws))
    ## } else {
    return(list(vs = vs, error = error, ks = ks))
    ##}
}
exact <- list(find_exact(1, y, Y, Ds), find_exact(2, y, Y, Ds))
## get_vs(exact)
## get_vs(solns[1:2])
all.equal(get_vs(exact), get_vs(solns[1:2]))

                                        # GBB
calc_obj(as.numeric(GBB), GBB.obs)
                                        # DART
calc_obj(as.numeric(DART), DART.obs)


                                        # Calculate ε for D <= 0.6
Ds <- localsolver::.doublewell$Ds
Dsr <- Ds[which(Ds < 0.6)]
idx <- seq_along(Dsr)
GBBr <- GBB[idx, ]
GBB.obsr <- GBB.obs[idx]
DARTr <- DART[idx, ]
DART.obsr <- DART.obs[idx]
Yr <- Y[idx, ]
yr <- y[idx]

calc_obj(Yr[, solns[[1]]$vs], yr)#, Yr)
sapply(2:4, function(i) calc_obj(rowMeans(Yr[, solns[[i]]$vs]), yr))
calc_obj(GBBr, GBB.obsr)#, Yr)
calc_obj(DARTr, DART.obsr)#, Yr)

                                        # What degree sequence does the n=4 node set have?
solns[[4]]$ks

rm(list = ls())

## What proportion of the dolphin n=4 opt sets have smaller error than $fixed and $rand?
ns <- readRDS("../data/ns-dolphin_doublewell.rds")
errors <- lapply(ns, get_error)
                                        # degree-preserving
sum(errors$fixed > max(errors$opt))
                                        # completely random
sum(errors$rand > max(errors$opt))

                                        # What degree do the BA network hub nodes have?
BA <- readRDS("../data/ba.rds")
sort(igraph::degree(BA), decreasing = TRUE)[1:5]

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

    return(list(vs = vs, error = error, ks = ks))
}
exact <- list(find_exact(1, y, Y, Ds), find_exact(2, y, Y, Ds))
all.equal(get_vs(exact), get_vs(solns[1:2]))

                                        # GBB
calc_obj(as.numeric(GBB), GBB.obs)
calc_obj(as.numeric(GBB), y)
                                        # DART
calc_obj(as.numeric(DART), DART.obs)
calc_obj(as.numeric(DART), y)


                                        # Calculate ε for D <= 0.6
Ds <- sdn::.doublewell$Ds
Dsr <- Ds[which(Ds < 0.6)]
idx <- seq_along(Dsr)
GBBr <- GBB[idx, ]
GBB.obsr <- GBB.obs[idx]
DARTr <- DART[idx, ]
DART.obsr <- DART.obs[idx]
Yr <- Y[idx, ]
yr <- y[idx]

calc_obj(Yr[, solns[[1]]$vs], yr)
sapply(2:4, function(i) calc_obj(rowMeans(Yr[, solns[[i]]$vs]), yr))
calc_obj(GBBr, GBB.obsr)
calc_obj(GBBr, yr)
calc_obj(DARTr, DART.obsr)
calc_obj(DARTr, yr)

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


                                        # Across all node set types (random, degree-preserving, and optimized) and
                                        # dynamics, what proportion of random node sets are larger than the largest
                                        # optimzied node set? and same for degree preserving.
count_smaller <- function(dynamics) {
    networks <- c(
        "dolphin", "celegans", "proximity", "euroroad", "email", "er", "ba", "hk", "gkk", "lfr",
        "drosophila", "reactome", "route_views", "spanish", "foldoc", "tree_of_life", "word_assoc",
        "enron", "marker_cafe", "prosper"        
    )
    allns <- lapply(
        networks,
        function(network) readRDS(paste0("../data/ns-", network, "_", dynamics, ".rds"))
    )
    allerrors <- lapply(allns, function(ns) lapply(ns, get_error))

    oe <- lapply(allerrors, `[[`, "opt")
    re <- lapply(allerrors, `[[`, "rand")
    fe <- lapply(allerrors, `[[`, "fixed")
    ##print(lengths(list(oe, re, fe)))

    c(
        rand = sum(mapply(function(opt, rand) sum(rand < max(opt)), oe, re)),
        fixed = sum(mapply(function(opt, fixed) sum(fixed < max(opt)), oe, fe))
    )
}

dyns <- c("doublewell", "mutualistic", "genereg", "SIS")
1 - rowSums(sapply(dyns, count_smaller)/2000)

#### Testing with DART and GBB
res <- data.frame(
    SNobs = y,
    DART = as.numeric(DART),
    DARTobs = DART.obs,
    GBB = as.numeric(GBB),
    GBBobs = GBB.obs
)

palette("Paired")
matplot(Ds, res[, c("SNobs", "DART", "DARTobs", "GBB", "GBBobs")],
        type = "l", lty = 1, col = c("black", 1:4), lwd = 5)
legend("topleft", bty = "n", lty = 1, col = c("black", 1:4), lwd = 5,
       legend = c("Unweighted mean state", "DART approx.", "DART obs.", "GBB approx.", "GBB obs."))

res$DARTint <- with(res, cumsum((DART - DARTobs)^2)/(length(DARTobs)*mean(DARTobs)))
res$DARTintSN <- with(res, cumsum((DART - SNobs)^2)/(length(SNobs)*mean(SNobs)))

matplot(Ds, res[, c("DARTint", "DARTintSN")], type = "l", col = c("red", "blue"), lty = 1, lwd = 5,
        xlab = "D", ylab = "Cumulative approximation error")
legend("bottomright", bty = "n", lty = 1, col = c("red", "blue"), lwd = 5,
       legend = c("DART against its own observable", "DART against the unweighted mean"))

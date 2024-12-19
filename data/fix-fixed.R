library(parallel)
ncores <- detectCores()-1
library(igraph)
library(optNS)

networks <- c(
    "dolphin", "celegans", "proximity", "euroroad", "email", "er", "gkk", "ba", "hk", "lfr",
    "drosophila", "reactome", "route_views", "spanish", "foldoc", "tree_of_life", "word_assoc", "enron",
    "marker_cafe", "prosper"
)
dynamics <- c(
    "doublewell", "mutualistic", "SIS", "genereg"
)

ntrials <- 100
for(network in networks) {
    ## network <- "email"
    g <- readRDS(paste0(network, ".rds")) # "../data/",  # this file is in ./data/
    N <- vcount(g)
    n <- floor(log(N))
    fullstatefile <- paste0("fullstate-", network, ".rds")

    for(dynamic in dynamics) {
        ## dynamic <- "doublewell"
        nsfile <- paste0("ns-", network, "_", dynamic, ".rds")

        Y <- readRDS(fullstatefile)[[dynamic]]
        y <- rowMeans(Y)

        ns <- readRDS(nsfile)
        best <- ns$opt[[which.min(get_error(ns$opt))]]

        fixeds <- make_dataset(
            ntrials = ntrials, ns.type = "fixed", ncores = ncores,
            n = n, g = g, comps = best$vs, y = y, Y = Y, optimize_weights = FALSE#optimize_weights
        )

        ns$fixed <- fixeds

        ## check
        test <- t(get_ks(ns$fixed))
        with(list(x = test[!duplicated(test), ]), {
            stopifnot(is.numeric(x) & length(x) == n)
        })

        saveRDS(ns, nsfile)
    }
}

## conds <- expand.grid(networks, dynamics, stringsAsFactors = FALSE)
## colnames(conds) <- c("network", "dynamics")

## nsfiles <- apply(conds, 1, function(row) {
##     row <- as.list(row)
##     paste0("ns-", row$network, "_", row$dynamics, ".rds")
## })

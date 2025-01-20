## This file adds a new column to the nodefeatures-*.csv files with the average error for each node for ~20 trials on each node set.

library(parallel)
ncores <- detectCores()-1
library(igraph)
library(optNS)

nets <- c(
    "dolphin", "proximity", "celegans", "euroroad", "email",
    "drosophila", "reactome", "route_views", "spanish", "foldoc",
    "tree_of_life", "word_assoc", "enron", "marker_cafe", "prosper",
    "er", "ba", "hk", "gkk", "lfr"
)
networks <- lapply(nets, function(network) {
    readRDS(paste0("../data/", network, ".rds"))
})
nodefeatures <- lapply(nets, function(network) {
    read.csv(paste0("../data/nodefeatures-", network, ".csv"))
})
fullstates <- lapply(nets, function(network) {
    readRDS(paste0("../data/fullstate-", network, ".rds"))
})
names(networks) <- names(nodefeatures) <- names(fullstates) <- nets

sample_size <- function(N, n) ceiling(20*N/n)

Ns <- sapply(networks, vcount)

for(net in nets) {
    print(net)
    g <- networks[[net]]
    nf <- nodefeatures[[net]]
    Y <- fullstates[[net]]$doublewell
    y <- rowMeans(Y)

    N <- Ns[net]
    n <- floor(log(N))
    samplesize <- sample_size(N, n)

    system.time(
        ns <- make_dataset(ntrials = samplesize, ns.type = "rand", ncores = ncores, n = n, g = g, y = y, Y = Y)
    )

    ns.vs <- t(get_vs(ns))
    dw.error <- get_error(ns)
    df <- data.frame(cbind(ns.vs, dw.error))

                                        # need to get the other four errors now
    df$mut.error <- with(list(Y = fullstates[[net]]$mutualistic), {
        y <- rowMeans(Y)
        apply(ns.vs, 1, function(vs) calc_obj(rowMeans(Y[, vs]), y))
    })
    df$SIS.error <- with(list(Y = fullstates[[net]]$SIS), {
        y <- rowMeans(Y)
        apply(ns.vs, 1, function(vs) calc_obj(rowMeans(Y[, vs]), y))
    })
    df$gene.error <- with(list(Y = fullstates[[net]]$genereg), {
        y <- rowMeans(Y)
        apply(ns.vs, 1, function(vs) calc_obj(rowMeans(Y[, vs]), y))
    })

                                        # reshape both the errors (to dynamics)
                                        # and nodes (to v, for the averaging)
    rdf1 <- reshape(
        df, varying = c("dw.error", "mut.error", "SIS.error", "gene.error"),
        v.names = "error",
        timevar = "dynamics",
        times = c("doublewell", "mutualistic", "SIS", "genereg"),
        new.row.names = seq(nrow(df)*4),
        direction = "long"
    )
    rdf2 <- reshape(
        rdf1, varying = paste0("V", seq(n)),
        v.names = "v", timevar = "selected", times = paste0("V", seq(n)),
        new.row.names = seq(nrow(rdf1)*n),
        direction = "long"
    )
                                        # and aggregate, taking the mean
    vdf <- aggregate(error ~ v + dynamics, data = rdf2, FUN = mean)
    
                                        # now, merge to nf
    dat <- merge(nf, vdf, by = c("v", "dynamics"), all.x = TRUE) # NB: this "error" is the mean error for each node

                                        # now, random forest
                                        # Have to save it so can use skl.
                                        # Going to add this new column to the CSV files that already exist
    write.csv(dat, paste0("../data/nodefeatures-", net, ".csv"), row.names = FALSE)
}

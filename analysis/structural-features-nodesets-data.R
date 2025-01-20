library(parallel)
ncores <- detectCores()-1
library(igraph)
library(optNS)

set.seed(123) # because of the Louvain algorithm
save_plots <- FALSE # TRUE
useall <- "no" # "yes"    # If no, only opt and rand
useweights <- "no" # "yes"
weightsflag <- switch(useweights, no = FALSE, yes = TRUE)
networks <- c(
    "dolphin", "proximity", "celegans", "euroroad", "email",
    "drosophila", "reactome", "route_views", "spanish", "foldoc",
    "tree_of_life", "word_assoc", "enron", "marker_cafe", "prosper",
    "er", "ba", "hk", "gkk", "lfr"
)

## If all cols in a df are the same class, apply()'s working object is a named vector of that class

get_nss <- function(net, dyn, ns.type) {
    nsname <- paste(net, dyn, sep = "_")
    nodesets[[nsname]][[ns.type]]
}

get_average_nodefeatures <- function(net, dyn, ns.type) {
    nss <- get_nss(net, dyn, ns.type) # the 100 node sets to work with
    nf <- nodefeatures[[net]]
    t(sapply(nss, function(ns) colMeans(nf[nf$v %in% ns$vs, -which(colnames(nf) == "v")])))
}

                                        # helper function
max_fromsame <- function(ns, part) {
    mbr <- part[ns$vs]
    max(table(mbr))
}

                                        # max number of nodes coming from the same community
get_pairs_scores <- function(net, dyn, ns.type) {
    nss <- get_nss(net, dyn, ns.type)
    mbr <- memberships[[net]]
    sapply(nss, max_fromsame, mbr)
}

                                        # average shortest path distance between nodes in a node set
get_avg_geods <- function(net, dyn, ns.type) { 
    vs <- t(get_vs(get_nss(net, dyn, ns.type)))
    apply(vs, 1, function(v) {
        geod <- distances(graphlist[[net]], v, v)
        mean(geod[lower.tri(geod)])
    })
}

graphlist <- lapply(networks, function(network) readRDS(paste0("../data/", network, ".rds")))
## fullstates <- lapply(networks, function(network) readRDS(paste0("../data/fullstate-", network, ".rds")))
nodefeatures <- lapply(networks, function(network) {
    nf <- read.csv(paste0("../data/nodefeatures-", network, ".csv"))
    cols <- c("v", "k", "cc", "bc", "knn", "lcl", "kcore")
    subset(nf, dynamics == "doublewell", select = cols)
})
names(graphlist) <- names(nodefeatures) <- networks # <- names(fullstates) 

dynamics <- c("doublewell", "SIS", "mutualistic", "genereg")

ns.types <- c("opt", "fixed", "rand", "constr", "quant", "knnconstr", "comm")
ns.types <- ns.types[switch(useall, no = c(1:3), yes = 1:length(ns.types))]            # change here

conds <- expand.grid(networks, dynamics, ns.types, stringsAsFactors = FALSE)
colnames(conds) <- c("net", "dyn", "ns.type")
nslistnames <- apply(conds[, 1:2], 1, paste, collapse = "_")
nodesets <- lapply(nslistnames, function(cond) {
    readRDS(paste0("../data/ns-", cond, switch(useweights, no = ".rds", yes = "_w.rds")))
}) # _w if useweights
names(nodesets) <- nslistnames

                                        # CHANGING VARIABLE NAMES
system.time(partitions <- lapply(graphlist, cluster_louvain)) # 37.840
memberships <- lapply(partitions, membership)
modularities <- sapply(partitions, modularity)
names(partitions) <- names(memberships) <- names(modularities) <- networks

## Separately, need to collect all of the node data from the nodefeatures csvs and get ready for the model
## k, bc, cc, knn, lcl, kcore; then community membership and average shortest path distance
system.time(dfs <- mcmapply(function(net, dyn, ns.type) { # 500 seconds
    nfs <- get_average_nodefeatures(net, dyn, ns.type)
    pairs <- get_pairs_scores(net, dyn, ns.type)
    geods <- get_avg_geods(net, dyn, ns.type)
    data.frame(
        network = net, dynamics = dyn, ns.type = ns.type, nfs, pairs = pairs, geods = geods,
        row.names = seq_len(nrow(nfs))
    )
}, conds$net, conds$dyn, conds$ns.type, SIMPLIFY = FALSE, mc.cores = ncores))
df <- do.call(rbind, dfs)
rownames(df) <- seq_len(nrow(df))

errors <- apply(
    conds, 1,
    function(row) get_error(get_nss(row["net"], row["dyn"], row["ns.type"])),
    simplify = FALSE
)
df$error <- do.call(c, errors)

write.csv(df, "../data/nodeset-features.csv", row.names = FALSE)

library(parallel)
ncores <- detectCores()-1
library(igraph)

useall <- "no" # "yes"
useweights <- "no" # "yes"
weightsflag <- switch(useweights, no = FALSE, yes = TRUE)

## Can I improve the get_* functions? I want to avoid global variables
get_from_ns <- function(row, ref, varname, collection) {
    network <- row[1]
    dynamic <- row[2]
    ns.type <- row[3]
    nslist <- paste(c(network, dynamic), collapse = "_")

    ns <- collection[[nslist]][[ns.type]]
    ref <- ref[[network]]

    var <- sapply(ns, function(x) mean(ref[x$vs], na.rm = TRUE))
    df <- data.frame(network = network, dynamics = dynamic, ns.type = ns.type, var = var,
                     row.names = NULL)
    colnames(df)[4] <- varname
    return(df)
}

count_pairs <- function(ns, part) {## nodeset, partition membership
    mbr <- part[ns$vs]
    A <- outer(mbr, mbr, `==`)
    h <- graph_from_adjacency_matrix(A, "undirected", diag = FALSE)
    ecount(h)
}

get_pairs_scores <- function(row, collection, partitions) {
    network <- row[1]
    dynamics <- row[2]
    ns.type <- row[3]
    nslist <- paste(c(network, dynamics), collapse = "_")
    partition <- partitions[[network]]

    ns <- collection[[nslist]][[ns.type]]

    var <- sapply(ns, count_pairs, partition)
    data.frame(network = network, dynamics = dynamics, ns.type = ns.type, pairs = var, row.names = NULL)
}

networks <- c( # only the empirical networks
    "dolphin", "celegans", "proximity", "euroroad", "email"
)
graphlist <- lapply(networks, function(network) readRDS(paste0("../data/", network, ".rds")))
dynamics <- c("doublewell", "SIS", "mutualistic", "genereg")
ns.types <- c("opt", "fixed", "rand", "constr", "quant", "comm")[switch(useall, no = 1:3, yes = 1:6)]
conds <- expand.grid(networks, dynamics, ns.types, stringsAsFactors = FALSE)
colnames(conds) <- c("networks", "dynamics", "ns.type")
nslistnames <- apply(conds[, 1:2], 1, paste, collapse = "_")

fullstates <- lapply(networks, function(network) readRDS(paste0("../data/fullstate-", network, ".rds")))
names(fullstates) <- networks
nodesets <- lapply(nslistnames, function(cond) {
    readRDS(paste0("../data/ns-", cond, switch(useweights, no = ".rds", yes = "_w.rds")))
}) # _w if useweights
names(nodesets) <- nslistnames

ks <- lapply(graphlist, degree)
knns <- lapply(graphlist, function(g) knn(g)$knn)
lcls <- lapply(graphlist, transitivity, type = "localundirected")
partitions <- lapply(graphlist, function(g) membership(cluster_louvain(g)))
names(ks) <- names(knns) <- names(lcls) <- names(partitions) <- networks


allk <- do.call(rbind, apply(conds, 1, get_from_ns, ks, "k", nodesets, simplify = FALSE))
allknn <- do.call(rbind, apply(conds, 1, get_from_ns, knns, "knn", nodesets, simplify = FALSE))
alllcl <- do.call(rbind, apply(conds, 1, get_from_ns, lcls, "lcl", nodesets, simplify = FALSE))
allpairs <- do.call(rbind, apply(conds, 1, get_pairs_scores, nodesets, partitions, simplify = FALSE))

df <- cbind(allk, allknn$knn, alllcl$lcl, allpairs$pairs)
colnames(df) <- c("network", "dynamics", "ns.type", "k", "knn", "lcl", "pairs")

df$network <- factor(df$network) # celegans is the reference. That's probably ok, as it approximates a BA network
df$dynamics <- factor(df$dynamics) # dw is the reference. That is ok.
df$ns.type <- factor(df$ns.type) # opt is the reference. I think I want rand to be the reference
df$ns.type <- relevel(df$ns.type, "rand") # make it so

## model is going to be
## mk <- glm(k ~ dynamics + network + ns.type, family = gaussian("log"), data = df)
## update(mk, . ~ ., + network:ns.type)
## each df needs to be network, dynamics, ns.type, [var], where var %in% k, knn, lcl, pairs

mk <- glm(# for flexibility
    k ~ dynamics + network + ns.type,
    family = gaussian("log"),
    data = df
)
mki <- update(mk, . ~ . + network:ns.type)

mknn <- glm(
    knn ~ dynamics + network + ns.type,
    family = gaussian("log"),
    data = df
)
mknni <- update(mknn, . ~ . + network:ns.type)

mlcl <- glm(
    sqrt(lcl) ~ dynamics + network + ns.type,
    family = gaussian,
    data = df
)
mlcli <- update(mlcl, . ~ . + network:ns.type)

mp <- glm(
    pairs ~ dynamics + network + ns.type,
    family = poisson,
    data = df
)
mpi <- update(mp, . ~ . + network:ns.type)

lapply(list(k_only = mk, knn = mknn, clustering = mlcl, communities = mp), summary)
lapply(list(k_only = mki, knn = mknni, clustering = mlcli, communities = mpi), summary)
    

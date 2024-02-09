library(parallel)
ncores <- detectCores()-1
library(igraph)

source("../src/functions.R")
get_error <- function(dl) sapply(dl, `[[`, "error")
load_fullstate <- function(net) {
    ## this is needed because all of the `fullstate` saved objects have the same name once loaded
    load(paste0("../data/fullstate-", net, ".rda"))
    temp <- fullstate
    rm(fullstate)
    return(temp)
}

empiricals <- c("dolphin", "celegans", "proximity", "euroroad", "email")
networks <- c(
    "dolphin", "celegans", "proximity", "euroroad", "email", "lfr"
)
for(net in networks) load(paste0("../data/", net, ".rda")); rm(net)
fullstates <- lapply(networks, load_fullstate)
names(fullstates) <- networks
dynamics <- c("dw", "SIS", "mutualistic", "genereg")

conds <- expand.grid(networks, dynamics)
colnames(conds) <- c("networks", "dynamics")

alloptnames <- apply(conds, 1, function(row) paste(row, collapse = "_"))
allopts <- mcmapply(
    function(net, dyn) readRDS(paste0("../data/optimized-nodesets/", paste(c(net, dyn), collapse = "-"), ".rds")),
    conds$networks, conds$dynamics, SIMPLIFY = FALSE, mc.cores = ncores
)
names(allopts) <- alloptnames

generate_nodesets <- function(network, dynamic) {
    source("../src/functions.R", local = TRUE)
    ## To do this I need the correct set of bparam values, y, Y, and the network
    g <- upgrade_graph(get(network))
    N <- vcount(g)
    k <- degree(g)
    Y <- fullstates[[network]][[dynamic]]
    y <- rowMeans(Y)
    bparam <- switch(
        dynamic,
        dw = doublewell_parms$Ds,
        SIS = SIS_parms$Ds,
        mutualistic = mutualistic_parms$Ds,
        genereg = genereg_parms$Ds
    )
    
    cond <- paste(c(network, dynamic), collapse = "_")
    opt <- allopts[[cond]]
    ntrials <- length(opt)

    bestopt <- opt[[which.min(get_error(opt))]]
    n <- length(bestopt$vs)

    fixed <- make_dataset(n, ntrials, bparam, y, Y, comps = bestopt$vs)
    rand <- make_dataset(n, ntrials, bparam, y, Y)
    constr <- make_dataset(n, ntrials, bparam, y, Y, use_connections = TRUE)
    quant <- make_dataset(n, ntrials, bparam, y, Y, use_connections = TRUE, use_quantiles = TRUE)

    return(list(opt = opt, fixed = fixed, rand = rand, constr = constr, quant = quant))
}

nodesets <- apply(conds, 1, function(row) generate_nodesets(row[1], row[2]), simplify = FALSE)
names(nodesets) <- alloptnames

get_partition <- function(net) {
    g <- upgrade_graph(get(net))
    communities <- cluster_fast_greedy(g)
    communities$membership
}
partitions <- lapply(networks, get_partition)
names(partitions) <- networks

count_pairs <- function(ns, part) {## nodeset, partition membership
    mbr <- part[ns$vs]
    A <- outer(mbr, mbr, `==`)
    h <- graph_from_adjacency_matrix(A, "undirected", diag = FALSE)
    ecount(h)
}

get_pairs_scores <- function(net, dyn) {
    cond <- paste(c(net, dyn), collapse = "_")
    nsc <- nodesets[[cond]]
    partition <- partitions[[net]]
    nsnames <- names(nsc)

    dl <- mapply(function(ns, nomen) {
        data.frame(pairs = sapply(ns, count_pairs, partition), network = net, dynamics = dyn, ns.type = nomen,
                   row.names = NULL)
    }, nsc, nsnames, SIMPLIFY = FALSE)

    df <- do.call(rbind, dl)
    rownames(df) <- seq_len(nrow(df))
    return(df)
}


## ok, so that's the data. Need to get it from each network.
allpairs <- apply(conds, 1, function(row) get_pairs_scores(row[1], row[2]), simplify = FALSE)

df <- do.call(rbind, allpairs)
df$network <- factor(df$network) # default ref is celegans
df$dynamics <- factor(df$dynamics) # default ref is dw
df$ns.type <- factor(df$ns.type)
df$ns.type <- relevel(df$ns.type, "opt")

model.aov <- aov(pairs ~ network + dynamics + ns.type, data = df)
model.glm <- glm(pairs ~ network + dynamics + ns.type, data = df, family = poisson)
exp(model.glm$coefficients)
summary(model.glm)
TukeyHSD(model.aov, "ns.type")

panelA <- hist(subset(df, network == "proximity" & dynamics == "dw" & ns.type == "rand")$pairs, breaks = 0:10, plot = FALSE)
panelB <- hist(subset(df, network == "proximity" & dynamics == "dw" & ns.type == "opt")$pairs, breaks = 0:10, plot = FALSE)
ylim <- range(unlist(c(panelA$counts, panelB$counts)))
dev.new(width = 14)
par(mfrow = c(1, 2))
plot(panelA, ylim = ylim, main = "Random", xlab = "Number of pairs from same community")
plot(panelB, ylim = ylim, main = "Optimized", xlab = "Number of pairs from same community")


test <- aggregate(pairs ~ network + dynamics + ns.type, FUN = mean, data = subset(df, ns.type %in% c("rand", "opt")))
reshape(test, timevar = "ns.type", idvar = c("network", "dynamics"), direction = "wide")

df$network <- relevel(df$network, "lfr")
model.glm2 <- glm(
    pairs ~ dynamics + network + ns.type,
    data = subset(df, ns.type %in% c("rand", "opt")),
    family = poisson
)
model.glm2int <- glm(
    pairs ~ dynamics + network * ns.type,
    data = subset(df, ns.type %in% c("rand", "opt")),
    family = poisson
)
## AIC(model.glm2, model.glm2int)
## summary(model.glm2)
summary(model.glm2int)

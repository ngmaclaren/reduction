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
    "dolphin", "celegans", "proximity", "euroroad", "email", "er", "gkk", "ba", "hk", "lfr"
)
for(net in networks) load(paste0("../data/", net, ".rda"))
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

get_lcl <- function(net) {
    g <- upgrade_graph(get(net))
    transitivity(g, "localundirected")
}
lcls <- lapply(networks, get_lcl)
names(lcls) <- networks

get_ns_lcl <- function(net, dyn) {
    cond <- paste(c(net, dyn), collapse = "_")
    nsc <- nodesets[[cond]]
    lcl <- lcls[[net]]
    nsnames <- names(nsc)

    dl <- mapply(function(ns, nomen) {
        data.frame(lcl = sapply(ns, function(S) mean(lcl[S$vs], na.rm = TRUE)),
                   network = net, dynamics = dyn, ns.type = nomen,
                   row.names = NULL)
    }, nsc, nsnames, SIMPLIFY = FALSE)

    df <- do.call(rbind, dl)
    rownames(df) <- seq_len(nrow(df))
    return(df)
}

alllcl <- apply(conds, 1, function(row) get_ns_lcl(row[1], row[2]), simplify = FALSE)

df <- do.call(rbind, alllcl)

df$network <- factor(df$network) # default ref is celegans
df$dynamics <- factor(df$dynamics) # default ref is dw
df$ns.type <- factor(df$ns.type)
df$ns.type <- relevel(df$ns.type, "opt")

hist(subset(df, network == "ba")$lcl)

model.glm <- glm(
    lcl ~ dynamics + network + ns.type, family = gaussian,
    data = subset(df, network %in% empiricals & ns.type %in% c("rand", "opt"))
)
model.glmint <- glm(
    lcl ~ dynamics + network * ns.type, family = gaussian,
    data = subset(df, network %in% empiricals & ns.type %in% c("rand", "opt"))
)
summary(model.glm)
summary(model.glmint)

model.aov <- aov(
    lcl ~ network + dynamics + ns.type,
    data = subset(df, network %in% empiricals)
)
TukeyHSD(model.aov, "ns.type")

hist(subset(df, network %in% empiricals)

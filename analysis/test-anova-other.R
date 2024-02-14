## This result needs to be on the /foreign/ dynamics.
## Set up should be similar, but there are more dynamics to compare.

library(parallel)
ncores <- detectCores()-1
library(igraph)

source("../src/functions.R")
get_error <- function(dl) sapply(dl, `[[`, "error")
load_fullstate <- function(net) {
    load(paste0("../data/fullstate-", net, ".rda"))
    temp <- fullstate
    rm(fullstate)
    return(temp)
}

useall <- "no"

empiricals <- c("dolphin", "celegans", "proximity", "euroroad", "email")
networks <- c(
    "dolphin", "celegans", "proximity", "euroroad", "email", "er", "gkk", "ba", "hk", "lfr"
)
for(net in networks) load(paste0("../data/", net, ".rda"))
fullstates <- lapply(networks, load_fullstate)
names(fullstates) <- networks
dynamics <- c("dw", "SIS", "mutualistic", "genereg")

                                        # conditions for testing foreign dynamics
conds <- expand.grid(networks, dynamics, dynamics)
colnames(conds) <- c("network", "dynamicsA", "dynamicsB")
conds <- conds[-which(conds$dynamicsA == conds$dynamicsB), ]

                                        # gather all of the opt node sets
optconds <- expand.grid(networks, dynamics)
colnames(optconds) <- c("networks", "dynamics")
alloptnames <- apply(optconds, 1, function(row) paste(row, collapse = "_"))
allopts <- mcmapply(
    function(net, dyn) readRDS(paste0("../data/optimized-nodesets/", paste(c(net, dyn), collapse = "-"), ".rds")),
    optconds$networks, optconds$dynamics, SIMPLIFY = FALSE, mc.cores = ncores
)
names(allopts) <- alloptnames

                                        # gather all of the random node sets, too
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
    if(useall == "yes") {
        constr <- make_dataset(n, ntrials, bparam, y, Y, use_connections = TRUE)
        quant <- make_dataset(n, ntrials, bparam, y, Y, use_connections = TRUE, use_quantiles = TRUE)
        comm <- make_dataset(n, ntrials, bparam, y, Y, use_connections = TRUE, use_communities = TRUE)

        return(list(opt = opt, fixed = fixed, rand = rand, constr = constr, quant = quant, comm = comm))
    } else {
        return(list(opt = opt, fixed = fixed, rand = rand))
    }
}

nodesets <- apply(optconds, 1, function(row) generate_nodesets(row[1], row[2]), simplify = FALSE)
names(nodesets) <- alloptnames

## now, we need to go through conds (with dynamicsA and dynamicsB) and evaluate each net/dynA combo against each other dynB
## use net and dynA to select the right node sets from `nodesets`, then evaluate against fullstate[[dynB]].
evaluate_errors <- function(conditions) {
    net <- conditions[1]#as.character(conditions$network)
    dynA <- conditions[2]#as.character(conditions$dynamicsA)
    dynB <- conditions[3]#as.character(conditions$dynamicsB)

    nsc <- nodesets[[paste(c(net, dynA), collapse = "_")]]
    nsnames <- names(nsc)
    
    g <- upgrade_graph(get(net))
    N <- vcount(g)
    k <- degree(g)
    Y <- fullstates[[net]][[dynB]]
    y <- rowMeans(Y)
    bparam <- switch(
        dynB,
        dw = doublewell_parms$Ds,
        SIS = SIS_parms$Ds,
        mutualistic = mutualistic_parms$Ds,
        genereg = genereg_parms$Ds
    )

    dl <- mapply(function(ns, nomen) {
        error <- sapply(ns, function(S) obj_fn(S$vs, y, Y, bparam))
        data.frame(
            error = error, network = net, dynamicsA = dynA, dynamicsB = dynB, ns.type = nomen, row.names = NULL
        )
    }, nsc, nsnames, SIMPLIFY = FALSE)

    df <- do.call(rbind, dl)
    rownames(df) <- seq_len(nrow(df))
    return(df)
}

## this is a list of foreign dynamics errors. It is the response variable.
alldynBerrors <- apply(conds, 1, evaluate_errors, simplify = FALSE)
## need to turn it into a df with the appropriate predictors: net, dynA, dynB, ns.type
df <- do.call(rbind, alldynBerrors)

df$network <- factor(df$network)
df$dynamicsA <- factor(df$dynamicsA)
df$dynamicsB <- factor(df$dynamicsB)
df$ns.type <- factor(df$ns.type)
df$ns.type <- relevel(df$ns.type, "rand")

model.glm <- glm(
    log(error) ~ network + dynamicsA + dynamicsB + ns.type,
    data = df,
    ## data = subset(df, network %in% empiricals),
    family = gaussian
)
summary(model.glm)
CI.glm <- confint(model.glm)
exp(CI.glm["(Intercept)", ])
exp(CI.glm["(Intercept)", ] + CI.glm["ns.typeopt", ])
exp(CI.glm["(Intercept)", ] + CI.glm["ns.typefixed", ])

model.aov <- aov(
    log(error) ~ network + dynamicsA + dynamicsB + ns.type,
    data = df
    ## data = subset(df, network %in% empiricals)
)
HSD <- TukeyHSD(model.aov, "ns.type")
show(HSD)

beta <- coefficients(model.glm)

## mean error for random node sets, across all networks and dynamics (A and B)
exp(beta[1])

## mean error for optimized node sets
exp(beta[1] + beta["ns.typeopt"])
exp(beta[1] + HSD$ns.type["opt-rand", "lwr"])
exp(beta[1] + HSD$ns.type["opt-rand", "upr"])

## mean error for fixed-degree node sets
exp(beta[1] + beta["ns.typefixed"])
exp(beta[1] + HSD$ns.type["fixed-rand", "lwr"])
exp(beta[1] + HSD$ns.type["fixed-rand", "upr"])

library(parallel)
ncores <- detectCores()-1

library(igraph)

source("../src/functions.R")

load_fullstate <- function(net) {
    ## this is needed because all of the `fullstate` saved objects have the same name once loaded
    load(paste0("../data/fullstate-", net, ".rda"))
    temp <- fullstate
    rm(fullstate)
    return(temp)
}

get_error <- function(dl) sapply(dl, `[[`, "error")
## get_vs? _ks?

## this version of generate_nodesets only needs rand, so maybe I can just select randomly without the big function
## No, I still want it, because it calculates the error
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
    n <- length(opt[[1]]$vs)

    rand <- make_dataset(n, ntrials, bparam, y, Y)

    return(list(opt = opt, rand = rand))
}

                                        # the "get_" functions
get_k <- function(net) {
    g <- upgrade_graph(get(net))
    degree(g)
}

get_ns_k <- function(net, dyn) {
    cond <- paste(c(net, dyn), collapse = "_")
    nsc <- nodesets[[cond]]
    k <- ks[[net]]
    nsnames <- names(nsc)

    dl <- mapply(function(ns, nomen) {
        data.frame(k = sapply(ns, function(S) mean(k[S$vs], na.rm = TRUE)), # mean vs. max
                   network = net, dynamics = dyn, ns.type = nomen,
                   row.names = NULL)
    }, nsc, nsnames, SIMPLIFY = FALSE)

    df <- do.call(rbind, dl)
    rownames(df) <- seq_len(nrow(df))
    return(df)
}

get_knn <- function(net) {
    g <- upgrade_graph(get(net))
    knn(g)$knn
}

get_ns_knn <- function(net, dyn) {
    cond <- paste(c(net, dyn), collapse = "_")
    nsc <- nodesets[[cond]]
    knn <- knns[[net]]
    nsnames <- names(nsc)

    dl <- mapply(function(ns, nomen) {
        data.frame(knn = sapply(ns, function(S) mean(knn[S$vs], na.rm = TRUE)),
                   network = net, dynamics = dyn, ns.type = nomen,
                   row.names = NULL)
    }, nsc, nsnames, SIMPLIFY = FALSE)

    df <- do.call(rbind, dl)
    rownames(df) <- seq_len(nrow(df))
    return(df)
}

get_lcl <- function(net) {
    g <- upgrade_graph(get(net))
    transitivity(g, "localundirected")
}

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

get_partition <- function(net) {
    g <- upgrade_graph(get(net))
    ##communities <- cluster_fast_greedy(g)
    communities <- cluster_louvain(g)
    membership(communities)
}

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


                                        # load the networks and other support vars
networks <- c( # only the empirical networks
    "dolphin", "celegans", "proximity", "euroroad", "email"
)
for(net in networks) load(paste0("../data/", net, ".rda"))
fullstates <- lapply(networks, load_fullstate)
names(fullstates) <- networks
dynamics <- c("dw", "SIS", "mutualistic", "genereg")

                                        # combinations of network and dynamics
conds <- expand.grid(networks, dynamics)
colnames(conds) <- c("networks", "dynamics")
alloptnames <- apply(conds, 1, function(row) paste(row, collapse = "_"))

                                        # load all optimized node sets
allopts <- mcmapply(
    function(net, dyn) readRDS(paste0("../data/optimized-nodesets/", paste(c(net, dyn), collapse = "-"), ".rds")),
    conds$networks, conds$dynamics, SIMPLIFY = FALSE, mc.cores = ncores
)
names(allopts) <- alloptnames

                                        # GLOBAL VARIABLES!!!
                                        # Make the nodesets list
nodesets <- apply(conds, 1, function(row) generate_nodesets(row[1], row[2]), simplify = FALSE)
names(nodesets) <- alloptnames
                                        # and the other needed variables
                                        # As a side note, this is why using global variables is bad practice.
ks <- lapply(networks, get_k)
names(ks) <- networks
knns <- lapply(networks, get_knn)
names(knns) <- networks
lcls <- lapply(networks, get_lcl)
names(lcls) <- networks
partitions <- lapply(networks, get_partition)
names(partitions) <- networks



## these are all lists of data frames, collect with do.call(rbind, ...)
allk <- do.call(rbind, apply(conds, 1, function(row) get_ns_k(row[1], row[2]), simplify = FALSE))
allknn <- do.call(rbind, apply(conds, 1, function(row) get_ns_knn(row[1], row[2]), simplify = FALSE))
alllcl <- do.call(rbind, apply(conds, 1, function(row) get_ns_lcl(row[1], row[2]), simplify = FALSE))
allpairs <- do.call(rbind, apply(conds, 1, function(row) get_pairs_scores(row[1], row[2]), simplify = FALSE))

## Merge not working. Because of the way I made these, I can use cbind.
df <- cbind(allk, allknn$knn, alllcl$lcl, allpairs$pairs)
colnames(df) <- c("k", "network", "dynamics", "ns.type", "knn", "lcl", "pairs") # fragile!
df <- df[, c("network", "dynamics", "ns.type", "k", "knn", "lcl", "pairs")]
df$network <- factor(df$network) # celegans is the reference. That's probably ok, as it approximates a BA network
df$dynamics <- factor(df$dynamics) # dw is the reference. That is ok.
df$ns.type <- factor(df$ns.type) # opt is the reference. I think I want rand to be the reference
df$ns.type <- relevel(df$ns.type, "rand")

## Now, the models.
## No need for TukeyHSD here because there is only one contrast of interest, between opt and rand.

## About a 10% difference overall, but network matters. 
mk <- glm(# for flexibility
    k ~ dynamics + network + ns.type,
    family = gaussian("log"),
    data = df
)
mki <- update(mk, . ~ . + network:ns.type)

                                        # 2/8/24: removed k as a predictor from all of the below
## makes a significant difference, but not a practical difference. About 2.5% difference. 
mknn <- glm(
    knn ~ dynamics + network + ns.type,
    family = gaussian("log"),
    data = df
)
mknni <- update(mknn, . ~ . + network:ns.type)

## no real effect on local clustering based on ns.type
mlcl <- glm(
    sqrt(lcl) ~ dynamics + network + ns.type,
    family = gaussian,
    data = df
)
mlcli <- update(mlcl, . ~ . + network:ns.type)

## only makes a difference for the dolphin network, which is by far the smallest network and has the fewest communities
mp <- glm(
    pairs ~ dynamics + network + ns.type,
    family = poisson,
    data = df
)
mpi <- update(mp, . ~ . + network:ns.type)

lapply(list(k_only = mk, knn = mknn, clustering = mlcl, communities = mp), summary)

lapply(list(k_only = mki, knn = mknni, clustering = mlcli, communities = mpi), summary)

### To work out the interaction:
## ns.type is a dummy variable. The reference category is rand, so the coefficient "ns.typeopt" expresses the change in the conditional mean (i.e., after accounting for the dyanmics and network dummies) between rand and opt.
## With the interaction involved it's more complicated. The base is now celegans and rand, when we go to a different network we add that network coefficient, e.g. networkdolphin. When we go to opt we add ns.typeopt. For different network and type opt we add (Intercept) + networkdolphin + ns.typeopt + networkdolphin:ns.typeot

## Let's try to look at some interaction plots.
## set dynamics==dw
themodel <- mlcli
thelines <- with(list(x = coefficients(themodel)), {
    beta0 <- x["(Intercept)"]
    beta1 <- x["ns.typeopt"]
    list(
        celegans = c(beta0, beta0 + beta1),
        dolphin = c(beta0 + x["networkdolphin"],
                    beta0 + x["networkdolphin"] + beta1 + x["networkdolphin:ns.typeopt"]),
        email = c(beta0 + x["networkemail"],
                  beta0 + x["networkemail"] + beta1 + x["networkemail:ns.typeopt"]),
        euroroad = c(beta0 + x["networkeuroroad"],
                  beta0 + x["networkeuroroad"] + beta1 + x["networkeuroroad:ns.typeopt"]),
        proximity = c(beta0 + x["networkproximity"],
                      beta0 + x["networkproximity"] + beta1 + x["networkproximity:ns.typeopt"])
    )
})
thesegments <- with(list(x = confint(themodel)), {
    beta0 <- x["(Intercept)", ]
    beta1 <- x["ns.typeopt", ]
    list(
        celegans = list(
            lower = c(beta0[1], beta0[1] + beta1[1]),
            upper = c(beta0[2], beta0[2] + beta1[2])
        ),
        dolphin = list(
            lower = c(beta0[1] + x["networkdolphin", 1],
                      beta0[1] + x["networkdolphin", 1] + beta1[1] + x["networkdolphin:ns.typeopt", 1]),
            upper = c(beta0[2] + x["networkdolphin", 2],
                      beta0[2] + x["networkdolphin", 2] + beta1[2] + x["networkdolphin:ns.typeopt", 2])
        ),
        email = list(
            lower = c(beta0[1] + x["networkemail", 1],
                      beta0[1] + x["networkemail", 1] + beta1[1] + x["networkemail:ns.typeopt", 1]),
            upper = c(beta0[2] + x["networkemail", 2],
                      beta0[2] + x["networkemail", 2] + beta1[2] + x["networkemail:ns.typeopt", 2])
        ),
        euroroad = list(
            lower = c(beta0[1] + x["networkeuroroad", 1],
                      beta0[1] + x["networkeuroroad", 1] + beta1[1] + x["networkeuroroad:ns.typeopt", 1]),
            upper = c(beta0[2] + x["networkeuroroad", 2],
                      beta0[2] + x["networkeuroroad", 2] + beta1[2] + x["networkeuroroad:ns.typeopt", 2])
        ),
        proximity = list(
            lower = c(beta0[1] + x["networkproximity", 1],
                      beta0[1] + x["networkproximity", 1] + beta1[1] + x["networkproximity:ns.typeopt", 1]),
            upper = c(beta0[2] + x["networkproximity", 2],
                      beta0[2] + x["networkproximity", 2] + beta1[2] + x["networkproximity:ns.typeopt", 2])
        )
    )
})
xlim <- c(1, 2)
ylim <- exp(range(unlist(c(thelines, thesegments))))

plot(NULL, xlim = xlim, ylim = ylim, xlab = "", xaxt = "no", ylab = "Mean")
for(i in seq_along(thelines)) {
    points(1:2, exp(thelines[[i]]), col = i, pch = 0)
    lines(1:2, exp(thelines[[i]]), col = i)
    segments(x0 = 1, y0 = exp(thesegments[[i]]$lower[1]), y1 = exp(thesegments[[i]]$upper[1]), col = i)
    segments(x0 = 2, y0 = exp(thesegments[[i]]$lower[2]), y1 = exp(thesegments[[i]]$upper[2]), col = i)
}
axis(1, at = 1:2, labels = c("Random", "Optimized"))
legend(
    "topleft", bty = "n", ncol = length(networks), col = seq_along(networks), pch = 0, legend = levels(df$network)
)



## celegans has a very large CV, and is the only one where the degree is substantially smaller in the optimized node sets.
##   dolphin  celegans proximity  euroroad     email 
## 0.5763020 1.8688823 0.6233392 0.4784803 0.9710590 
## I wonder if I can find something similar for pairs. Maybe when there is a very large modularity?
## The dolphin network doesn't have a particularly high modularity
##   dolphin  celegans proximity  euroroad     email 
## 0.5201139 0.4324883 0.7097282 0.8672712 0.5717054
## It does have relatively few communities
## dolphin  celegans proximity  euroroad     email 
##       4        11         6        22        11 
## And it is a small network.

## I think the quantiling thing may be a red herring. But does it matter? Should we quantile without removing the very largest nodes (except in extremely heterogeneous networks?)

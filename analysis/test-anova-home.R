## This result is on the /home/ dynamics. I still need to do the /foreign/ dynamics.
## What I'm seeing here is that opt is always better and it's better than everything. Of course.
## Random is always worst.
## Fixed-degree, constrained, and quantiled are statistically the same as each other. This means we don't on average gain anything by taking the extra quantiling step. However, we do gain something on average by rejecting the largest degree nodes. I need to check this again after I've fixed the make_dataset() function. 

## if(interactive()) {
##     if(getwd() != "/user/neilmacl/Documents/reduction/analysis/") {
##         setwd("/user/neilmacl/Documents/reduction/analysis/")
##     }
## }

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
    ## quant <- make_dataset(n, ntrials, bparam, y, Y, use_connections = TRUE, use_quantiles = TRUE)
    comm <- make_dataset(n, ntrials, bparam, y, Y, use_connections = TRUE, use_communities = TRUE)

    return(list(opt = opt, fixed = fixed, rand = rand, constr = constr, comm = comm)) # quant = quant
}

nodesets <- apply(conds, 1, function(row) generate_nodesets(row[1], row[2]), simplify = FALSE)
names(nodesets) <- alloptnames

collect_errors <- function(network, dynamic) { # node set collection; list of opt, fixed, rand, constr, quant
    ## I want to return a data frame with error, network, dynamics, node set
    cond <- paste(c(network, dynamic), collapse = "_")
    nsc <- nodesets[[cond]]
    nsnames <- names(nsc)

    dl <- mapply(function(ns, nomen) {
        error <- get_error(ns)
        data.frame(error = error, network = network, dynamics = dynamic, ns.type = nomen, row.names = NULL)
    }, nsc, nsnames, SIMPLIFY = FALSE)

    df <- do.call(rbind, dl)
    rownames(df) <- seq_len(nrow(df))
    return(df)
}

allerrors <- apply(conds, 1, function(row) collect_errors(row[1], row[2]), simplify = FALSE)

df <- do.call(rbind, allerrors)
df$network <- factor(df$network)
df$dynamics <- factor(df$dynamics)
df$ns.type <- factor(df$ns.type)
df$ns.type <- relevel(df$ns.type, "opt")

model.glm <- glm(
    log(error) ~ network + dynamics + ns.type,
    data = df,
    family = gaussian
)
summary(model.glm)

                                        # this is a different model b/c log(y) != log(mean(y) + ε)
model.aov <- aov(log(error) ~ network + dynamics + ns.type, data = df)
TukeyHSD(model.aov, "ns.type") # the p-value is already adjusted for multiple comparisons. How? Doesn't say in the docs.


## panel it by something. Dynamics?
palette("Set 1")
dev.new(width = 16, height = 8)
par(mfrow = c(2, 2))
for(i in seq_along(levels(df$dynamics))) {
    plotdf <- subset(df, dynamics == levels(df$dynamics)[i])
    means <- aggregate(error ~ network + ns.type, FUN = mean, data = plotdf)
    plot(
        x = jitter(as.numeric(plotdf$network), amount = 0.1), xlab = "Network",
        y = plotdf$error, ylab = "Error",
        col = adjustcolor(as.numeric(plotdf$ns.type), 0.5), pch = 16, cex = 0.5,
        axes = FALSE, ylim = c(min(means$error)*.9, max(means$error)*1.05)#log = "y"
    )
    points(
        x = as.numeric(means$network),
        y = means$error,
        col = as.numeric(means$ns.type), bg = adjustcolor(as.numeric(means$ns.type), .5), pch = 21, cex = 3
    )
    box()
    axis(1, labels = levels(plotdf$network), at = seq_along(levels(plotdf$network)))
    axis(2)
    legend("bottomright", bty = "n", col = seq_along(levels(plotdf$ns.type)), pch = 16,
           legend = levels(plotdf$ns.type))
}

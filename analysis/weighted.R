library(parallel)
ncores <- detectCores()-1

## Load the same random data sets from before, as well as the opts.
networks <- c(
    "dolphin", "celegans", "proximity", "euroroad", "email", "er", "gkk", "ba", "hk", "lfr"
)
dynamics <- c("doublewell", "SIS", "mutualistic", "genereg")
ns.types <- c("opt", "rand")
conds <- expand.grid(networks, dynamics, ns.types, stringsAsFactors = FALSE)
colnames(conds) <- c("networks", "dynamics", "ns.type")
nslistnames <- apply(conds[!duplicated(conds[, 1:2]), 1:2], 1, paste, collapse = "_")
fullstates <- lapply(networks, function(network) readRDS(paste0("../data/fullstate-", network, ".rds")))
names(fullstates) <- networks
unweighted_nodesets <- lapply(nslistnames, function(cond) {
    readRDS(paste0("../data/ns-", cond, ".rds"))
}) # _w if useweights
names(unweighted_nodesets) <- nslistnames

## then, load the opts from the weighted case
weighted_nodesets <- lapply(nslistnames, function(cond) readRDS(paste0("../data/ns-", cond, "_w.rds")))
names(weighted_nodesets) <- nslistnames

## need to do both the home analysis and the other analysis
collect_errors <- function(row, collection) { # this will be apply()ed to the `conds` df
    ## output: data.frame(error = error, network = network, dynamics = dynamics, ns.type = ns.type)
    network <- row[1]
    dynamic <- row[2]
    ns.type <- row[3]
    nslist <- paste(c(network, dynamic), collapse = "_")
    
    data.frame(
        error = optNS::get_error(collection[[nslist]][[ns.type]]),
        network = network,
        dynamics = dynamic,
        ns.type = ns.type,
        row.names = NULL
    )
}

udf <- do.call(rbind, apply(conds, 1, collect_errors, collection = unweighted_nodesets))
wdf <- do.call(rbind, apply(conds, 1, collect_errors, collection = weighted_nodesets))

wdf <- wdf[wdf$ns.type == "opt", ]
wdf$ns.type[wdf$ns.type == "opt"] <- "optw"

df <- rbind(udf, wdf)

df$network <- factor(df$network)
df$dynamics <- factor(df$dynamics)
df$ns.type <- factor(df$ns.type)
df$ns.type <- relevel(df$ns.type, "rand")

model.glm <- glm(
    log(error) ~ network + dynamics + ns.type,
    data = subset(df, network != "er"),
    family = gaussian
)
summary(model.glm)

## Ok, that's pretty much it there.
## Now for the other dynamics.
conds_other <- expand.grid(networks, dynamics, dynamics, ns.types, stringsAsFactors = FALSE)
colnames(conds_other) <- c("networks", "dynamicsA", "dynamicsB", "ns.type")

collect_errors_other <- function(row, collection, weightflag) {
    network <- row[1]
    dynamicsA <- row[2]
    dynamicsB <- row[3]
    ns.type <- row[4]
    nslist <- paste(c(network, dynamicsA), collapse = "_")

    Y <- fullstates[[network]][[dynamicsB]]
    y <- rowMeans(Y)

    error <- sapply(
        collection[[nslist]][[ns.type]],
        function(x) optNS::obj_fn(x$vs, y = y, Y = Y,
                                  optimize_weights = switch(ns.type, rand = FALSE, opt = weightflag),
                                  ws = x$ws)
    )

    data.frame(
        error = error, network = network, dynamicsA = dynamicsA, dynamicsB = dynamicsB, ns.type = ns.type,
        row.names = NULL
    )
}

udfo <- do.call(
    rbind,
    mclapply(
        split(conds_other, seq(nrow(conds_other))),
        function(row) collect_errors_other(as.character(row), unweighted_nodesets, FALSE),
        mc.cores = ncores
    )
)

wdfo <- do.call(
    rbind,
    mclapply(
        split(conds_other, seq(nrow(conds_other))),
        function(row) collect_errors_other(as.character(row), weighted_nodesets, TRUE),
        mc.cores = ncores
    )
)

wdfo <- wdfo[wdfo$ns.type == "opt", ]
wdfo$ns.type[wdfo$ns.type == "opt"] <- "optw"

dfo <- rbind(udfo, wdfo)
dups <- apply(dfo, 1, function(row) any(duplicated(row)))
dfo <- dfo[!dups, ]

dfo$network <- factor(dfo$network)
dfo$dynamicsA <- factor(dfo$dynamicsA)
dfo$dynamicsB <- factor(dfo$dynamicsB)
dfo$ns.type <- factor(dfo$ns.type)
dfo$ns.type <- relevel(dfo$ns.type, "rand")

model.glm <- glm(
    log(error) ~ network + dynamicsA + dynamicsB + ns.type,
    data = subset(dfo, network != "er"),
    family = gaussian
)
summary(model.glm)

1/exp(coef(model.glm)[c("ns.typeopt", "ns.typeoptw")])

model.aov <- aov(log(error) ~ network + dynamicsA + dynamicsB + ns.type, data = subset(dfo, network != "er"))
anova(model.aov)
TukeyHSD(model.aov, "ns.type")
## drat. Need to go back and confirm now, though: didn't catch before that I was including self-dynamics in the other-dynamics analysis.

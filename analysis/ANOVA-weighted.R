library(parallel)
ncores <- detectCores()-1

networks <- c(
                                        # Exclude ER
    "dolphin", "proximity", "celegans", "euroroad", "email", "gkk", "ba", "hk", "lfr" # , "er"
)
                                        # Use all dynamics
dynamics <- c("doublewell", "SIS", "mutualistic", "genereg")
                                        # Want to compare unweighted random, unweighted optimized,
                                        # and weight-optimized
ns.types <- c("opt", "rand")

names.ns.all <- apply(expand.grid(networks, dynamics), 1, function(row) paste(row, collapse = "_"))

                                        # Need all full states to recompute errors
fullstates <- lapply(networks, function(network) readRDS(paste0("../data/fullstate-", network, ".rds")))
names(fullstates) <- networks

                                        # Unweighted node sets for random and opt
nodesets <- lapply(names.ns.all, function(nsname) readRDS(paste0("../data/ns-", nsname, ".rds")))
                                        # For weight-optimized node sets
nodesets_w <- lapply(names.ns.all, function(nsname) readRDS(paste0("../data/ns-", nsname, "_w.rds")))
names(nodesets) <- names(nodesets_w) <- names.ns.all

### Home dynamics ###

conds <- expand.grid(networks, dynamics, ns.types, stringsAsFactors = FALSE)
colnames(conds) <- c("network", "dynamics", "ns.type")

collect_errors <- function(row, nodesets, fullstates) {
    network <- row[1]
    dynamic <- row[2]
    ns.type <- row[3]

    if(ns.type == "opt") {
        nsname <- paste(c(network, dynamic), collapse = "_")
                                        # Uses the stored error for opt ns, including ns-*_w (which was computed with optimize_weights = TRUE)
        error <- optNS::get_error(nodesets[[nsname]][[ns.type]])
    } else if(ns.type == "rand") {
        nsname <- paste(c(network, "doublewell"), collapse = "_")

        Y <- fullstates[[network]][[dynamic]]
        y <- rowMeans(Y)

        ns <- nodesets[[nsname]][[ns.type]]
                                        # Only gets called if rand, so optimze_weights = FALSE is correct
                                        # Will actually do this for the ns-*_w rand sets, too, but we will discard those
        error <- sapply(ns, function(x) optNS::obj_fn(x$vs, y = y, Y = Y, optimize_weights = FALSE))
    }

    data.frame(
        error = error, network = network, dynamics = dynamic, ns.type = ns.type,
        row.names = NULL
    )
}

                                        # u for unweighted
udf <- do.call(rbind, apply(conds, 1, collect_errors, nodesets, fullstates))
                                        # w for weighted; use the correct *_w node sets
wdf <- do.call(rbind, apply(conds, 1, collect_errors, nodesets_w, fullstates))
                                        # drop the rand ns
wdf <- wdf[wdf$ns.type == "opt", ]
                                        # and relable the opt ns
wdf$ns.type[wdf$ns.type == "opt"] <- "optw"

df <- rbind(udf, wdf)

df$network <- factor(df$network, levels = networks)
df$dynamics <- factor(df$dynamics, levels = dynamics)
df$ns.type <- factor(df$ns.type, levels = c("rand", "opt", "optw"))

## model.glm <- glm(log(error) ~ network + dynamics + ns.type, data = df, family = gaussian)
model.aov <- aov(log(error) ~ network + dynamics + ns.type, data = df)
model.lm <- lm(log(error) ~ network + dynamics + ns.type, data = df)
summary(model.lm)$r.squared
anova(model.aov)
## summary(model.glm)
## coef(model.glm)
## exp(coef(model.glm))
## 1/exp(coef(model.glm))
hsd <- TukeyHSD(model.aov)$ns.type
hsd
1/exp(hsd[, "diff"])

### Other dyanmics ###

conds_other <- expand.grid(networks, dynamics, dynamics, ns.types, stringsAsFactors = FALSE)
colnames(conds_other) <- c("networks", "dynamicsA", "dynamicsB", "ns.type")
dups <- apply(conds_other, 1, function(row) any(duplicated(row)))
conds_other <- conds_other[!dups, ]

collect_errors_other <- function(row, nodesets, fullstates, optimize_weights = FALSE) {
    network <- row[1]
    dynamicA <- row[2]
    dynamicB <- row[3]
    ns.type <- row[4]

                                        # error will be handled the same for all, so just need the right node set
    if(ns.type == "opt") {
                                        # use the node set for dynamics A
        nsname <- paste(c(network, dynamicA), collapse = "_")
    } else if(ns.type == "rand") {
                                        # for heur/rand, use the doublewell node sets
                                        # Dynamics doesn't matter for these, so use the reference/standard
        nsname <- paste(c(network, "doublewell"), collapse = "_")
    }

                                        # Evaluate error against dynamics B (should be no case where A=B)
    Y <- fullstates[[network]][[dynamicB]]
    y <- rowMeans(Y)

    ns <- nodesets[[nsname]][[ns.type]]

    error <- sapply(ns, function(x) {
                                        # Need to include the possibility that we'll want to use the stored ws, but not the stored error (because cross dynamics)
        optNS::obj_fn(x$vs, y = y, Y = Y, optimize_weights = optimize_weights, ws = x$ws)
    })

    data.frame(
        error = error, network = network, dynamicsA = dynamicA, dynamicsB = dynamicB, ns.type = ns.type,
        row.names = NULL
    )
}

udfo <- do.call(
    rbind,
    mclapply(
        split(conds_other, seq(nrow(conds_other))),
                                        # default is optimize_weights = FALSE
        function(row) collect_errors_other(as.character(row), nodesets, fullstates),
        mc.cores = ncores
    )
)

wdfo <- do.call(
    rbind,
    mclapply(
        split(conds_other, seq(nrow(conds_other))),
                                        # Use the stored weights to optimize
                                        # and make sure to use the nodesets_w, the ones with stored weights
        function(row) collect_errors_other(as.character(row), nodesets_w, fullstates, TRUE), # FALSE
        mc.cores = ncores
    )
)

wdfo <- wdfo[wdfo$ns.type == "opt", ]
wdfo$ns.type[wdfo$ns.type == "opt"] <- "optw"

dfo <- rbind(udfo, wdfo)

dfo$network <- factor(dfo$network, levels = networks)
dfo$dynamicsA <- factor(dfo$dynamicsA, levels = dynamics)
dfo$dynamicsB <- factor(dfo$dynamicsB, levels = dynamics)
dfo$ns.type <- factor(dfo$ns.type, levels = c("rand", "opt", "optw"))
dfo$ns.type <- relevel(dfo$ns.type, "rand")

##model.glm_other <- glm(log(error) ~ network + dynamicsA + dynamicsB + ns.type, data = dfo, family = gaussian)
model.aov_other <- aov(log(error) ~ network + dynamicsA + dynamicsB + ns.type, data = subset(dfo, network != "er"))
model.lm_other <- lm(log(error) ~ network + dynamicsA + dynamicsB + ns.type, data = subset(dfo, network != "er"))

##summary(model.glm_other)

##1/exp(coef(model.glm_other))

summary(model.lm_other)$r.squared
anova(model.aov_other)
TukeyHSD(model.aov_other, "ns.type")

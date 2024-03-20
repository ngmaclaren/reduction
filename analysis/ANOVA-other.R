                                        # Opens the same way as ANOVA-home, but the conditions will be different
networks <- c(
                                        # Exclude ER
    "dolphin", "celegans", "proximity", "euroroad", "email", "gkk", "ba", "hk", "lfr" # , "er"
)
                                        # Use all dynamics
dynamics <- c("doublewell", "SIS", "mutualistic", "genereg")
ns.types <- c("opt", "fixed", "rand", "constr", "quant", "knnconstr", "comm")
                                        # The different ns.types need different handling
                                        # These two we need for each dynamics
names.ns.opt <- c("opt", "fixed")
                                        # These we'll just use the [network]_doublewell node sets
names.ns.heur <- c("rand", "constr", "quant", "knnconstr", "comm")

names.ns.all <- apply(expand.grid(networks, dynamics), 1, function(row) paste(row, collapse = "_"))

                                        # Need all full states to recompute errors
fullstates <- lapply(networks, function(network) readRDS(paste0("../data/fullstate-", network, ".rds")))
names(fullstates) <- networks

nodesets <- lapply(names.ns.all, function(nsname) readRDS(paste0("../data/ns-", nsname, ".rds")))
names(nodesets) <- names.ns.all

                                        # Different from here
conds <- expand.grid(networks, dynamics, dynamics, ns.types, stringsAsFactors = FALSE)
colnames(conds) <- c("network", "dynamicsA", "dynamicsB", "ns.type")
                                        # exclude cases in which A=B
dups <- apply(conds, 1, function(row) any(duplicated(row)))
conds <- conds[!dups, ]

collect_errors <- function(row, nodesets, fullstates) {
    network <- row[1]
    dynamicA <- row[2]
    dynamicB <- row[3]
    ns.type <- row[4]

                                        # error will be handled the same for all, so just need the right node set
    if(ns.type %in% c("opt", "fixed")) {
                                        # use the node set for dynamics A
        nsname <- paste(c(network, dynamicA), collapse = "_")
    } else if(ns.type %in% c("rand", "constr", "quant", "knnconstr", "comm")) {
                                        # for heur/rand, use the doublewell node sets
                                        # Dynamics doesn't matter for these, so use the reference/standard
        nsname <- paste(c(network, "doublewell"), collapse = "_")
    }

                                        # Evaluate error against dynamics B (should be no case where A=B)
    Y <- fullstates[[network]][[dynamicB]]
    y <- rowMeans(Y)

    ns <- nodesets[[nsname]][[ns.type]]

    error <- sapply(ns, function(x) optNS::obj_fn(x$vs, y = y, Y = Y, optimize_weights = FALSE))

    data.frame(
        error = error, network = network, dynamicsA = dynamicA, dynamicsB = dynamicB, ns.type = ns.type,
        row.names = NULL
    )
}

df <- do.call(
    rbind,
    apply(conds, 1, function(row) collect_errors(as.character(row), nodesets, fullstates))
)

df$network <- factor(df$network)
df$dynamicsA <- factor(df$dynamicsA)
df$dynamicsB <- factor(df$dynamicsB)
df$ns.type <- factor(df$ns.type)
df$ns.type <- relevel(df$ns.type, "rand")

model.glm <- glm(log(error) ~ network + dynamicsA + dynamicsB + ns.type, data = df, family = gaussian)
model.aov <- aov(log(error) ~ network + dynamicsA + dynamicsB + ns.type, data = df)

anova(model.aov)
summary(model.glm)
TukeyHSD(model.aov, "ns.type")
1 - (model.glm$deviance/model.glm$null.deviance)
exp(coef(model.glm))
1/exp(coef(model.glm))


networks <- c(
                                        # Exclude ER
    "dolphin", "celegans", "proximity", "euroroad", "email", "gkk", "ba", "hk", "lfr" # , "er"
)
                                        # Use all dynamics
dynamics <- c("doublewell", "SIS", "mutualistic", "genereg")
ns.types <- c("opt", "fixed", "rand", "constr", "quant", "knnconstr", "comm")

names.ns.all <- apply(expand.grid(networks, dynamics), 1, function(row) paste(row, collapse = "_"))

                                        # Need all full states to recompute errors
fullstates <- lapply(networks, function(network) readRDS(paste0("../data/fullstate-", network, ".rds")))
names(fullstates) <- networks

nodesets <- lapply(names.ns.all, function(nsname) readRDS(paste0("../data/ns-", nsname, ".rds")))
names(nodesets) <- names.ns.all

conds <- expand.grid(networks, dynamics, ns.types, stringsAsFactors = FALSE)
colnames(conds) <- c("network", "dynamics", "ns.type")

collect_errors <- function(row, nodesets, fullstates) {
    network <- row[1]
    dynamic <- row[2]
    ns.type <- row[3]

                                        # this if else statement should output error
    if(ns.type %in% c("opt", "fixed")) {
        nsname <- paste(c(network, dynamic), collapse = "_")
        error <- optNS::get_error(nodesets[[nsname]][[ns.type]])
    } else if(ns.type %in% c("rand", "constr", "quant", "knnconstr", "comm")) {
                                        # For heuristic algs and completely random, use the doubewell node sets
                                        # This is arbitrary: node sets exist for each condition
        nsname <- paste(c(network, "doublewell"), collapse = "_")
                                        # But the error will need to be recomputed
        Y <- fullstates[[network]][[dynamic]]
        y <- rowMeans(Y)

        ns <- nodesets[[nsname]][[ns.type]]

        error <- sapply(ns, function(x) optNS::obj_fn(x$vs, y = y, Y = Y, optimize_weights = FALSE))
    }

    data.frame(
        error = error, network = network, dynamics = dynamic, ns.type = ns.type,
        row.names = NULL
    )
}

df <- do.call(
    rbind,
    apply(conds, 1, function(row) collect_errors(as.character(row), nodesets, fullstates))
)

df$network <- factor(df$network)
df$dynamics <- factor(df$dynamics)
df$ns.type <- factor(df$ns.type)
df$ns.type <- relevel(df$ns.type, "rand")

model.glm <- glm(log(error) ~ network + dynamics + ns.type, data = df, family = gaussian)
model.aov <- aov(log(error) ~ network + dynamics + ns.type, data = df)

anova(model.aov)
summary(model.glm)
TukeyHSD(model.aov, "ns.type")
1 - (model.glm$deviance/model.glm$null.deviance)
exp(coef(model.glm))
1/exp(coef(model.glm))

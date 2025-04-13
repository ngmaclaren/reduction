use_weighted_networks <- FALSE # TRUE
use_directed_networks <- FALSE # TRUE
exclude_heuristic <- TRUE # FALSE

networks <- c(
                                        # Exclude ER
    "dolphin", "proximity", "celegans", "euroroad", "email", "gkk", "ba", "hk", "lfr", # , "er"
    "drosophila", "reactome", "route_views", "spanish", "foldoc", "tree_of_life", "word_assoc",
    "enron", "marker_cafe", "prosper"
)
weighted <- c(
    "windsurfers", "macaques", "train_terrorists", "highschool", "drug", "residence_hall", "netsci_weighted",
    "proximity_weighted", "gap_junction_herm", "intl_trade"
)
directed <- c(
    "canton", "physician_trust", "email_company", "flamingo", "ecoli", "yeast", "usair", "jung-c", "email_uni",
    "faa"
)
if(use_weighted_networks) networks <- weighted
if(use_directed_networks) networks <- directed
                                        # Use all dynamics and node set types
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

if(exclude_heuristic) {
    retain <- c("opt", "fixed", "rand")
    conds <- conds[conds$ns.type %in% retain, ]
}

collect_errors <- function(row, nodesets, fullstates) {
    network <- row[1]
    dynamic <- row[2]
    ns.type <- row[3]

                                        # this if else statement should output approximation error
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

df$network <- factor(df$network, levels = networks)
df$dynamics <- factor(df$dynamics, levels = dynamics)
df$ns.type <- factor(df$ns.type, levels = c("rand", "opt", "fixed", "constr", "quant", "knnconstr", "comm"))
df$logerror <- log(df$error)

## model.aov <- aov(log(error) ~ network + dynamics + ns.type, data = df)
system.time(model.aov <- aov(logerror ~ network + dynamics + ns.type, data = df))
system.time(model.lm <- lm(logerror ~ network + dynamics + ns.type, data = df))

summary(model.lm)$r.squared
anova(model.aov)
## summary(model.glm)
round(TukeyHSD(model.aov, "ns.type")$ns.type, 3)
## 1 - (model.glm$deviance/model.glm$null.deviance)
## exp(coef(model.aov))
## 1/exp(coef(model.aov))
    
## exp(coef(model.glm))
## 1/exp(coef(model.glm))

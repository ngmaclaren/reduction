useall <- "no" # "yes"
useweights <- "no" # "yes"

networks <- c(
    "dolphin", "celegans", "proximity", "euroroad", "email", "er", "gkk", "ba", "hk", "lfr"
)
dynamics <- c("doublewell", "SIS", "mutualistic", "genereg")
ns.types <- c("opt", "fixed", "rand", "constr", "quant", "comm")[switch(useall, no = 1:3, yes = 1:6)]

conds <- expand.grid(networks, dynamics, ns.types, stringsAsFactors = FALSE)
colnames(conds) <- c("networks", "dynamics", "ns.type")

nslistnames <- apply(conds[!duplicated(conds[, 1:2]), 1:2], 1, paste, collapse = "_")

fullstates <- lapply(networks, function(network) readRDS(paste0("../data/fullstate-", network, ".rds")))
names(fullstates) <- networks
nodesets <- lapply(nslistnames, function(cond) {
    readRDS(paste0("../data/ns-", cond, switch(useweights, no = ".rds", yes = "_w.rds")))
}) # _w if useweights
names(nodesets) <- nslistnames

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

df <- do.call(rbind, apply(conds, 1, collect_errors, collection = nodesets))
df$network <- factor(df$network)
df$dynamics <- factor(df$dynamics)
df$ns.type <- factor(df$ns.type)
df$ns.type <- relevel(df$ns.type, "rand")

model.glm <- glm(
    log(error) ~ network + dynamics + ns.type,
    data = df,
    family = gaussian
)
summary(model.glm)

model.aov <- aov(log(error) ~ network + dynamics + ns.type, data = df)
TukeyHSD(model.aov, "ns.type")

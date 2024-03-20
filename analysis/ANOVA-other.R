library(parallel)
ncores <- detectCores()-1

useall <- "yes" # "no"
useweights <- "no" # "yes"
weightsflag <- switch(useweights, no = FALSE, yes = TRUE)

networks <- c(
    "dolphin", "celegans", "proximity", "euroroad", "email", "er", "gkk", "ba", "hk", "lfr"
)
dynamics <- c("doublewell", "SIS", "mutualistic", "genereg")
ns.types <- c("opt", "fixed", "rand", "constr", "quant", "knnconstr", "comm")
ns.types <- ns.types[switch(useall, no = 1:3, yes = 1:length(ns.types))]

conds <- expand.grid(networks, dynamics, dynamics, ns.types, stringsAsFactors = FALSE)
colnames(conds) <- c("networks", "dynamicsA", "dynamicsB", "ns.type")
dups <- apply(conds, 1, function(row) any(duplicated(row)))
conds <- conds[!dups, ]

## nslistnames <- apply(conds[!duplicated(conds[, 1:2]), 1:2], 1, paste, collapse = "_")
## nslistnames <- unique(apply(conds[, c("networks", "dynamicsA")], 1, paste, collapse = "_"))
nslistnames <- apply(expand.grid(networks, dynamics), 1, paste, collapse = "_")

fullstates <- lapply(networks, function(network) readRDS(paste0("../data/fullstate-", network, ".rds")))
names(fullstates) <- networks
nodesets <- lapply(nslistnames, function(cond) {
    readRDS(paste0("../data/ns-", cond, switch(useweights, no = ".rds", yes = "_w.rds")))
}) # _w if useweights
names(nodesets) <- nslistnames

collect_errors <- function(row, collection) {
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
                                  optimize_weights = switch(ns.type, FALSE, opt = weightsflag),
                                  ws = x$ws)
    )

    data.frame(
        error = error, network = network, dynamicsA = dynamicsA, dynamicsB = dynamicsB, ns.type = ns.type,
        row.names = NULL
    )
}

##df <- do.call(rbind, apply(conds, 1, collect_errors, collection = nodesets))
df <- do.call(
    rbind,
    mclapply(
        split(conds, seq(nrow(conds))),
        function(row) collect_errors(as.character(row), nodesets),
        ##collect_errors, collection = nodesets,
        mc.cores = ncores
    )
)
df$network <- factor(df$network)
df$dynamicsA <- factor(df$dynamicsA)
df$dynamicsB <- factor(df$dynamicsB)
df$ns.type <- factor(df$ns.type)
df$ns.type <- relevel(df$ns.type, "rand")

dups <- apply(df, 1, function(row) any(duplicated(row)))
df <- df[!dups, ]

model.glm <- glm(
    log(error) ~ network + dynamicsA + dynamicsB + ns.type,
    data = subset(df, network != "er"),
    family = gaussian
)
summary(model.glm)

model.aov <- aov(log(error) ~ network + dynamicsA + dynamicsB + ns.type, data = subset(df, network != "er"))
anova(model.aov)
TukeyHSD(model.aov, "ns.type")

## McFadden's pseudo-R2
1 - (model.glm$deviance/model.glm$null.deviance)

## pdiff <- function(m) {
##     x <- coefficients(m)
##     b <- x
##     b[-1] <- x[1] + x[-1]
##     b <- exp(b)
##     ((b[-1]/b[1])*100) - 100
## }

## ## sapply(modellist, pdiff)
## as.data.frame(pdiff(model.glm))
1/exp(coef(model.glm)["ns.typeopt"])
                                        # 14.8734 -- much more sane.
                                        # 7.094317 with weights. Wow: no reason at all to optimize weights.
                                        # Never mind: 38.39367 if I don't optimize the node weights of the random
                                        # node sets. Crazy. 

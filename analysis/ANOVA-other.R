## Based on ANOVA-home, but more complicated because need to evaluate on the other dynamics
library(parallel)
ncores <- detectCores()-1

useall <- "yes" # "no"
useweights <- "no" # "yes"
weightsflag <- switch(useweights, no = FALSE, yes = TRUE)

networks <- c(
    "dolphin", "celegans", "proximity", "euroroad", "email", "er", "gkk", "ba", "hk", "lfr"
)
dynamics <- c("doublewell", "SIS", "mutualistic", "genereg")
ns.types <- c("opt", "fixed", "rand", "constr", "quant", "comm")[switch(useall, no = 1:3, yes = 1:6)]

conds <- expand.grid(networks, dynamics, dynamics, ns.types, stringsAsFactors = FALSE)
colnames(conds) <- c("networks", "dynamicsA", "dynamicsB", "ns.type")

nslistnames <- apply(conds[!duplicated(conds[, 1:2]), 1:2], 1, paste, collapse = "_")

fullstates <- lapply(networks, function(network) readRDS(paste0("../data/fullstate-", network, ".rds")))
names(fullstates) <- networks
nodesets <- lapply(nslistnames, function(cond) {
    readRDS(paste0("../data/ns-", cond, switch(useweights, no = ".rds", yes = "_w.rds")))
}) # _w if useweights
names(nodesets) <- nslistnames

## need to evaluate errors before can collect them.
## alternately, maybe can do the evaluation inside?

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
        function(x) optNS::obj_fn(x$vs, y = y, Y = Y, optimize_weights = weightsflag, ws = x$ws)
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

model.glm <- glm(
    log(error) ~ network + dynamicsA + dynamicsB + ns.type,
    data = df,
    family = gaussian
)
summary(model.glm)

model.aov <- aov(log(error) ~ network + dynamicsA + dynamicsB + ns.type, data = df)
TukeyHSD(model.aov, "ns.type")


pdiff <- function(m) {
    x <- coefficients(m)
    b <- x
    b[-1] <- x[1] + x[-1]
    b <- exp(b)
    ((b[-1]/b[1])*100) - 100
}

## sapply(modellist, pdiff)
as.data.frame(pdiff(model.glm))

useall <- "yes" # "no"
useweights <- "no" # "yes"

networks <- c(
    "dolphin", "celegans", "proximity", "euroroad", "email", "er", "gkk", "ba", "hk", "lfr"
)
dynamics <- c("doublewell", "SIS", "mutualistic", "genereg")
ns.types <- c("opt", "fixed", "rand", "constr", "quant", "knnconstr", "comm")
ns.types <- ns.types[switch(useall, no = 1:3, yes = 1:length(ns.types))]

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
    data = subset(df, network != "er"),
    family = gaussian
)
summary(model.glm)

## model.lm <- lm(
##     log(error) ~ network + dynamics + ns.type,
##     data = subset(df, network != "er")
## )
## anova(model.lm)

model.aov <- aov(log(error) ~ network + dynamics + ns.type, data = subset(df, network != "er"))
anova(model.aov)
TukeyHSD(model.aov, "ns.type")

## McFadden's pseudo-R2
1 - (model.glm$deviance/model.glm$null.deviance)


## Below here is code to help interpret the coefficients of this model, and investigate alternatives.
## The simplest interpretation, though, is as follows:
## The coefficients on the ns.type levels express the average (across networks and dynamics) difference---on the log scale---between random node sets and the focal node set type. Improvement with optimization is dramatic:
1/exp(coef(model.glm)["ns.typeopt"])
                                        # 283.4513
                                        # 220.2669 # but we're optimizing node weights for random as well. Is that right?
## That's how many times smaller optimized node set approximation error is, on average.
1/exp(coef(model.glm)["ns.typefixed"])






##                                         # avg for opt
## b.opt <- exp(coefficients(model.glm)["(Intercept)"] + coefficients(model.glm)["ns.typeopt"])
##                                         # avg for degree-preserving
## b.fixed <- exp(coefficients(model.glm)["(Intercept)"] + coefficients(model.glm)["ns.typefixed"])
##                                         # avg for random
## b.rand <- exp(coefficients(model.glm)["(Intercept)"])

## 1 - (b.rand - b.fixed)/b.rand
## 1 - (b.rand - b.opt)/b.rand

## pdiff <- function(m) {
##     x <- coefficients(m)
##     b <- x
##     b[-1] <- x[1] + x[-1]
##     b <- exp(b)
##     ((b[-1]/b[1])*100) - 100
## }

## ## sapply(modellist, pdiff)
## as.data.frame(pdiff(model.glm))


## b <- coef(model.glm)
## exp(b["ns.typeopt"])
## 1/exp(b["ns.typeopt"]) # the mean optimized error is this many times smaller
## ## Verify?
## sapply(df[, 2:4], levels)
## test <- subset(df, network == "ba" & dynamics == "doublewell")
## exp(b["ns.typeopt"])*exp(mean(log(test$error[test$ns.type == "rand"])))
## mean(test$error[test$ns.type == "opt"])

## lapply(split(test, test$ns.type), function(x) summary(log(x$error)))
## 6.838-4.562

## summary(
##     glm(log(error) ~ ns.type, data = subset(df, network == "ba" & dynamics == "doublewell"), family = gaussian)
## )

## ## Ok, these agree, actually. The model coefficient of -5.647 must be the /average/ difference, on the log scale, between random and optimized node sets. Not all conditions achieve that difference, but we expect it to hold on average across the various contrasts.

## ## let's try a couple more
## lapply(networks[-6], function(net) {
##     summary(
##         glm(log(error) ~ ns.type,
##             data = subset(df, network == net & dynamics == "genereg"),
##             family = gaussian)
##     )
## })

## ## not sure what to do about this:
## summary(update(model.glm, . ~ . + dynamics:ns.type))
## ## seems to show that the actual error we achieve depends on the dynamics, but the overall pattern of results doesn't change.
## mean(c(-2.38, -2.38 + c(-6.149, -1.808, -5.098)))
##                                         # [1] -5.64375
## b["ns.typeopt"]
##                                         # -5.64704 # very close
## interactionmodel <- update(model.glm, . ~ . + dynamics:ns.type)
## sapply(list(model.glm, interactionmodel), function(m) 1 - (m$deviance/m$null.deviance))
##                                         # [1] 0.7012492 0.7955257
##                                         # explain 9.4% more variance by including the interaction term
## ## AIC(model.glm, interactionmodel)
## ##                  df      AIC
## ## model.glm        19 95466.31
## ## interactionmodel 37 85947.30

## ## Strangely, optimizing on the double-well model has the least improvement of any of the dynamics. Optimizing on genereg leads to dramatic improvement, for example. You can see these effects by looking at SI figs. S1--S3.

## ## Just because I'm curious...
## doubleinteraction <- update(interactionmodel, . ~ . + network:ns.type)
## ## relatively less improvement
## AIC(model.glm, interactionmodel, doubleinteraction)
## sapply(list(model.glm, interactionmodel, doubleinteraction), function(m) 1 - (m$deviance/m$null.deviance))
## ## So, knowing the dynamics is important. That definitely helps us know how much improvement we get by optimizing. However, the basic conclusion---that we get tremendous improvement by optimizing compared to any of the heuristic algorithms---is unaffected by the interaction. 

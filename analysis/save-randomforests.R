library(igraph)
library(randomForest)
library(optNS)

networks <- c(
    "dolphin", "proximity", "celegans", "euroroad", "email", "gkk", "ba", "hk", "lfr", "er",
    "drosophila", "reactome", "route_views",
    "spanish",
    "foldoc"
    ## "tree_of_life"
    ## "word_assoc"
    ## "enron", "marker_cafe", "prosper"
)
Ns <- sapply(networks, function(net) vcount(readRDS(paste0("../data/", net, ".rds"))))
Ms <- sapply(networks, function(net) ecount(readRDS(paste0("../data/", net, ".rds"))))

rfdf <- data.frame(network = networks, N = Ns, M = Ms)

for(network in networks) {
    network <- "drosophila"
    dynamics <- c("doublewell", "mutualistic", "SIS", "genereg")
    vsfile <- paste0("../data/nodefeatures-", network, ".rds")
    nsfiles <- paste0("../data/ns-", network, "_", dynamics, ".rds")

    vs <- readRDS(vsfile) # vertices
    ## ns <- readRDS(nsfile) # node sets
    ## ns <- do.call(c, lapply(nsfiles, function(nsfile) readRDS(nsfile)$opt))
    ns <- do.call(
        rbind,
        lapply(dynamics, function(dyn) {
            opts <- readRDS(paste0("../data/ns-", network, "_", dyn, ".rds"))$opt
            count <- as.data.frame(table(as.numeric(get_vs(opts))))
            colnames(count)[1] <- "v"
            count$dynamics <- dyn
            count
        })
    )

    ## count <- as.data.frame(table(as.numeric(get_vs(ns$opt))))
    ## count <- as.data.frame(table(as.numeric(get_vs(ns))))
    ## colnames(count)[1] <- "v"
    ## count$dynamics <- rep(dynamics, each = 100)

    ## df <- merge(vs, count, by = "v", all.x = TRUE)
    df <- merge(vs, ns, by = "v")
    df$Freq[is.na(df$Freq)] <- 0
    ## df$Sel <- as.numeric(df$Freq > 0)
    ## df$Freq <- factor(df$Freq, ordered = TRUE)
    ## df$Sel <- factor(df$Sel)

    tm <- system.time({
        rf <- randomForest(
            x = df[, c("k", "knn", "lcl", "cc", "bc", "kcore", "dynamics")], y = df[, "Freq"],
            ##Freq ~ k + knn + lcl + cc + bc + kcore,  data = df,
            importance = TRUE, proximity = FALSE, keep.forest = FALSE
        )
        ## if(length(unique(df$Sel)) > 1) { # doesn't work with dolphin because all nodes selected at least once
        ##     rf$sel <- randomForest(
        ##         x = df[, c("k", "knn", "lcl", "cc", "bc", "kcore")], y = df[, "Sel"],
        ##         ##Sel ~ k + knn + lcl + cc + bc + kcore,  data = df,
        ##         importance = TRUE, proximity = FALSE, keep.forest = FALSE
        ##     )
        ## } else {
        ##     rf$sel <- NULL
        ## }
    })

    rf <- rfPoisson(
        x = df[, c("k", "knn", "lcl", "cc", "bc", "kcore", "dynamics")], y = df[, "Freq"],
        importance = TRUE, proximity = FALSE, keep.forest = FALSE
    )
    
        

    saveRDS(rf.Sel, paste0("../data/randomforest-", network, ".rds"))
    
    cat(network, ", walltime = ", tm[3], "\n", sep = "")
}

## rfdf$walltime <- c(0.051, 0.597, 0.69, 2.202, 2.367, 1.834, 1.799, 1.737, 2.196, 2.238, rep(NA, 10))
## rfdf["drosophila", "walltime"] <- 8.887
## rfdf["reactome", "walltime"] <- 34.276
## rfdf["route_views", "walltime"] <- 40.507
## rfdf["spanish", "walltime"] <- 153.204
## rfdf["foldoc", "walltime"] <- 185.31
## rfdf["tree_of_life", "walltime"] <- 319.702
## rfdf["word_assoc", "walltime"] <- 820 # approximate
## lm(log(walltime) ~ log(N), data = rfdf)
## rfdf$predicted <- with(rfdf, exp(-10.165 + 1.648*log(N)))


### From here down, move to an analysis file

## the confusion matrix is at rf.Sel$confusion
## the rows are the actual condition, the columns are the predicted condition.
## In a binary classification, if we are trying to classify "0", the matrix is
## TP FN
## FP TN,
## where TP: true positive, FN: false negative, FP: false positive, TN: true negative
## Precision: TP/(TP + FP)
## Recall: TP/(TP + FN)
## F1 = (2TP)/(2TP + FP + FN)

show(rf$freq$confusion)
show(rf$freq$err.rate[nrow(rf$freq$err.rate), "OOB"])

cat(
    "Network: ", network, ", RF run time: ", tm[3], "\n\n", sep = ""
)

with(list(m = rf$freq$confusion), {
                                        # All based on out-of-bag error...
    TP <- m["1", "1"] # we selected "1" (column) and it was actually "1" (row)
    FP <- m["0", "1"] # we selected "1" (column) and it was actually "0" (row)
    FN <- m["1", "0"] # we selected "0" (column) and it was actually "1" (row)
    cat("Precision:", TP/(TP + FP), "\n",
        "Recall:", TP/(TP + FN), "\n",
        "F1:", 2*TP/(2*TP + FP + FN), "\n")
})

show(importance(rf$freq, type = 2))

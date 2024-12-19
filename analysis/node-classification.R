library(igraph)
library(randomForest)
library(optNS)

networks <- c(
  ##   "dolphin", "proximity", "celegans", "euroroad", "email", "gkk", "ba", "hk", "lfr", "er"
  ## , "drosophila", "reactome", "route_views"
  ## , "spanish", "foldoc"
  ## , "tree_of_life"
    ## "word_assoc"
    ## "enron"
    ## "marker_cafe"
    "prosper"
)
## Ns <- sapply(networks, function(net) vcount(readRDS(paste0("../data/", net, ".rds"))))
## Ms <- sapply(networks, function(net) ecount(readRDS(paste0("../data/", net, ".rds"))))

network <- "spanish"
dynamics <- "doublewell"
vsfile <- paste0("../data/nodefeatures-", network, ".rds")
nsfile <- paste0("../data/ns-", network, "_", dynamics, ".rds")

vs <- readRDS(vsfile)
ns <- readRDS(nsfile)

count <- as.data.frame(table(as.numeric(get_vs(ns$opt))))

df <- merge(vs, count, by.x = "v", by.y = "Var1", all.x = TRUE)
df$Freq[is.na(df$Freq)] <- 0
df$Freq <- factor(df$Freq, ordered = TRUE)
df$Sel <- as.numeric(df$Freq > 0)
df$Sel <- factor(df$Sel)

tm <- system.time(
    rf.Sel <- randomForest(Sel ~ k + knn + lcl + cc + bc + kcore,  data = df, importance = TRUE, proximity = TRUE)
)
## the confusion matrix is at rf.Sel$confusion
## the rows are the actual condition, the columns are the predicted condition.
## In a binary classification, if we are trying to classify "0", the matrix is
## TP FN
## FP TN,
## where TP: true positive, FN: false negative, FP: false positive, TN: true negative
## Precision: TP/(TP + FP)
## Recall: TP/(TP + FN)
## F1 = (2TP)/(2TP + FP + FN)

                                        # Only concerned with non-zero
## show(rf.Sel$err.rate[nrow(rf.Sel$err.rate), "OOB"])

cat(
    "Network: ", network, ", RF run time: ", tm[3], "\n\n", sep = ""
)

with(list(m = rf.Sel$confusion), {
                                        # All based on out-of-bag error...
    TP <- m["1", "1"] # we selected "1" (column) and it was actually "1" (row)
    FP <- m["0", "1"] # we selected "1" (column) and it was actually "0" (row)
    FN <- m["1", "0"] # we selected "0" (column) and it was actually "1" (row)
    cat("Precision:", TP/(TP + FP), "\n",
        "Recall:", TP/(TP + FN), "\n",
        "F1:", 2*TP/(2*TP + FP + FN), "\n")
})

show(importance(rf.Sel, type = 2))

library(igraph)

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

## nm <- data.frame(g = networks, N = Ns, M = Ms)

for (network in networks) {
    netfile <- paste0("../data/", network, ".rds")

    g <- readRDS(netfile)
    N <- vcount(g)
    M <- ecount(g)

    tm <- system.time(
        nodefeatures <<- data.frame( #! global assignment---shouldn't really be necessary but maybe it is
            v = as.numeric(V(g)),
            k = degree(g),
            knn = knn(g)$knn,
            lcl = transitivity(g, "localundirected", isolates = "zero"),
            cc = closeness(g, normalized = TRUE),
            bc = betweenness(g, normalized = TRUE),
            kcore = coreness(g)
        )
    )

    saveRDS(nodefeatures, paste0("../data/nodefeatures-", network, ".rds"))

    cat(network, ", N = ", N, ", M = ", M, ", walltime = ", tm[3], "\n", sep = "")
}

## nm$walltime <- c(0.004, 0.046, 0.045, 0.150, 0.320, 0.195, 0.184, 0.158, 0.184, 0.548, 1.553, 19.087, 7.977, 31.316, 74.914, 217.597, rep(NA, 4))
## lm(log(nm$walltime) ~ log(I(nm$N^3 * nm$M))) 
## nm$predicted <- exp(-14.57 + 0.88*log(nm$N*nm$M))
## nm$predicted <- exp(-15.1618 + 0.49*log(nm$N^3*nm$M))

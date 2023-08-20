library(igraph)

networks <- c(
    "dolphin", "proximity", "celegans", "euroroad", "email",
    "er", "ba", "hk", "gkk", "lfr"
)

for(net in networks) load(paste0("../data/", net, ".rda"))

report <- function(g) {
    k <- degree(g)
    cat(
        "N = ", vcount(g), "\n",
        "M = ", ecount(g), "\n",
        "k.bar = ", mean(k), "\n",
        "CV = ", mean(k)/sd(k), "\n",
        sep = ""
    )
}

for(net in networks) {
    print(net)
    report(get(net))
}

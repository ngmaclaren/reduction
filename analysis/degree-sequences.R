### BROKEN---DO NOT RUN


library(parallel)
ncores <- detectCores()-1
library(igraph)
library(optNS)

g <- readRDS("../data/dolphin.rds")
N <- vcount(g)
k <- degree(g)

Y <- readRDS("../data/fullstate-dolphin.rds")$doublewell
y <- rowMeans(Y)

ns <- 1:4
ntrials <- 100

opts <- lapply(ns, function(n) {
    make_dataset(
        ntrials = ntrials, ns.type = "opt", ncores = ncores,
        n = n, g = g, y = y, Y = Y
    )
})

poss <- seq_len(N)

lapply(opts, get_ks)



## test some ggplot histograms
x <- as.data.frame(t(lapply(opts, get_ks)[[4]]))
colnames(x) <- c("first", "second", "third", "fourth")
y <- reshape(x, varying = colnames(x), v.names = "k", timevar = "order", direction = "long")
y$order <- factor(y$order, ordered = TRUE)

library(ggplot2)
library(ggthemes)

dev.new(
ggplot(y, aes(k, fill = order)) +
    geom_histogram(binwidth = 1) +
    scale_fill_tableau(guide = "none") +
    theme_classic() +
    labs(x = "Degree", y = "Frequency")

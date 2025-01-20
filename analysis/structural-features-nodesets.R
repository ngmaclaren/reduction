## Add: average distance between nodes in a node set
## In terms of plots, I think all we're doing is community membership (K) and average distance

library(parallel)
ncores <- detectCores()-1
library(igraph)
library(optNS)

set.seed(123) # because of the Louvain algorithm
save_plots <- FALSE # TRUE
useall <- "no" # "yes"    # If no, only opt and rand
useweights <- "no" # "yes"
weightsflag <- switch(useweights, no = FALSE, yes = TRUE)
networks <- c(
    "dolphin", "proximity", "celegans", "euroroad", "email",
    "drosophila", "reactome", "route_views", "spanish", "foldoc",
    "tree_of_life", "word_assoc", "enron", "marker_cafe", "prosper",
    "er", "ba", "hk", "gkk", "lfr"
)

## If all cols in a df are the same class, apply()'s working object is a named vector of that class

get_nss <- function(net, dyn, ns.type) {
    nsname <- paste(net, dyn, sep = "_")
    nodesets[[nsname]][[ns.type]]
}

get_average_nodefeatures <- function(net, dyn, ns.type) {
    nss <- get_nss(net, dyn, ns.type) # the 100 node sets to work with
    nf <- nodefeatures[[net]]
    t(sapply(nss, function(ns) colMeans(nf[nf$v %in% ns$vs, -which(colnames(nf) == "v")])))
}

                                        # helper function
max_fromsame <- function(ns, part) {
    mbr <- part[ns$vs]
    max(table(mbr))
}

                                        # max number of nodes coming from the same community
get_pairs_scores <- function(net, dyn, ns.type) {
    nss <- get_nss(net, dyn, ns.type)
    mbr <- memberships[[net]]
    sapply(nss, max_fromsame, mbr)
}

                                        # average shortest path distance between nodes in a node set
get_avg_geods <- function(net, dyn, ns.type) { 
    vs <- t(get_vs(get_nss(net, dyn, ns.type)))
    apply(vs, 1, function(v) {
        geod <- distances(graphlist[[net]], v, v)
        mean(geod[lower.tri(geod)])
    })
}

graphlist <- lapply(networks, function(network) readRDS(paste0("../data/", network, ".rds")))
## fullstates <- lapply(networks, function(network) readRDS(paste0("../data/fullstate-", network, ".rds")))
nodefeatures <- lapply(networks, function(network) {
    nf <- read.csv(paste0("../data/nodefeatures-", network, ".csv"))
    cols <- c("v", "k", "cc", "bc", "knn", "lcl", "kcore")
    subset(nf, dynamics == "doublewell", select = cols)
})
names(graphlist) <- names(nodefeatures) <- networks # <- names(fullstates) 

dynamics <- c("doublewell", "SIS", "mutualistic", "genereg")

ns.types <- c("opt", "fixed", "rand", "constr", "quant", "knnconstr", "comm")
ns.types <- ns.types[switch(useall, no = c(1, 3), yes = 1:length(ns.types))]

conds <- expand.grid(networks, dynamics, ns.types, stringsAsFactors = FALSE)
colnames(conds) <- c("net", "dyn", "ns.type")
nslistnames <- apply(conds[, 1:2], 1, paste, collapse = "_")
nodesets <- lapply(nslistnames, function(cond) {
    readRDS(paste0("../data/ns-", cond, switch(useweights, no = ".rds", yes = "_w.rds")))
}) # _w if useweights
names(nodesets) <- nslistnames

                                        # CHANGING VARIABLE NAMES
system.time(partitions <- lapply(graphlist, cluster_louvain)) # 37.840
memberships <- lapply(partitions, membership)
modularities <- sapply(partitions, modularity)
names(partitions) <- names(memberships) <- names(modularities) <- networks

## Separately, need to collect all of the node data from the nodefeatures csvs and get ready for the model
## k, bc, cc, knn, lcl, kcore; then community membership and average shortest path distance
df <- apply(conds, 1, function(row) {
    net <- row["net"]
    dyn <- row["dyn"]
    ns.type <- row["ns.type"]
    nfs <- get_average_nodefeatures(net, dyn, ns.type)
    pairs <- get_pairs_scores(net, dyn, ns.type)
    geods <- get_avg_geods(net, dyn, ns.type)
    data.frame(
        network = net, dynamics = dyn, ns.type = ns.type, nfs, pairs = pairs, geods = geods
    )
})

## df is going to start with conds, then add the mean features for each 
## Use the histogram (more continuous) plot function below for avg shortest path distance
## Then use the more_discrete function for community membership
##


with(
    list(
        x = get_avg_geods("dolphin", "doublewell", "rand"),
        y = get_avg_geods("dolphin", "doublewell", "opt")
    ), {
        ## plot(x, y)
        hist.all <- hist(c(x, y), plot = FALSE)
        hist.x <- hist(x, breaks = hist.all$breaks, plot = FALSE)
        hist.y <- hist(y, breaks = hist.all$breaks, plot = FALSE)
        ylim <- range(c(hist.x$counts, hist.y$counts))
        plot(hist.x, col = adjustcolor(1, .5), ylim = ylim, xlab = "Average distance", main = "")
        plot(hist.y, col = adjustcolor(2, .5), ylim = ylim, add = TRUE)
    }
)






## 2025-01-07: I don't want to do these anymore. I want to load the csvs and get from there. 
ks <- lapply(graphlist, degree)
knns <- lapply(graphlist, function(g) knn(g)$knn)
lcls <- lapply(graphlist, transitivity, type = "localundirected")
partitions <- lapply(graphlist, function(g) membership(cluster_louvain(g)))
names(ks) <- names(knns) <- names(lcls) <- names(partitions) <- networks


allk <- do.call(rbind, apply(conds, 1, get_from_ns, ks, "k", nodesets, simplify = FALSE))
allknn <- do.call(rbind, apply(conds, 1, get_from_ns, knns, "knn", nodesets, simplify = FALSE))
alllcl <- do.call(rbind, apply(conds, 1, get_from_ns, lcls, "lcl", nodesets, simplify = FALSE))
allpairs <- do.call(rbind, apply(conds, 1, get_pairs_scores, nodesets, partitions, simplify = FALSE))

df <- cbind(allk, allknn$knn, alllcl$lcl, allpairs$pairs)
colnames(df) <- c("network", "dynamics", "ns.type", "k", "knn", "lcl", "pairs")


df$network <- factor(df$network, levels = networks)
df$dynamics <- factor(df$dynamics, levels = dynamics)
df$ns.type <- factor(df$ns.type, levels = c("rand", "opt", "fixed", "constr", "quant", "knnconstr", "comm"))

#### NEW MODEL # 2025-01-07: keeping the model, only plotting community and distance
moo <- glm(# I don't know why I came up with this name. 
    ns.type ~ dynamics + network + k + knn + lcl + pairs,
    family = binomial,
    data = subset(df, network != "er")
)
summary(moo)

                                        # uncomment to test dynamics interactions. Nothing to get excited about.
## summary(update(moo, . ~ . + dynamics:pairs))
print(1 - (moo$deviance/moo$null.deviance))

exp(coefficients(moo))
## So, odds ratio for knn is 1.02. For an increase of 1 in the knn, the odds of that node set being an optimized node set increase by a factor of 1.02.
                                        # percent change in odds (for an increase of one unit of x)
(exp(coefficients(moo)) - 1)*100
## These are all based on going from random to optimized. So, if we go from 0 pairs of nodes from the same community to 1 pair of nodes from the same community, it is 15% less likely that that node set is an optimized node set. Part of the reason for these differences in scale of the odds ratios, though, that there are large differences in scale in the predictor variables. Let's normalize the predictors. 

                                        # Can get to r2=0.03 or so. Not worth it.
                                        # results seem consistent, though
complicated <- glm(
    ns.type ~ dynamics + network + k + knn + lcl + pairs +
        network:k + network:knn + network:lcl + network:pairs,
    family = binomial, data = df
)
summary(complicated)
1 - (complicated$deviance/complicated$null.deviance)
  

## Plots here, at least for now, so I don't need ten different files. (That's an exaggeration, but still.)
## Density plots for knn and lcl. Barplot for pairs
## Select a network and a dynamics

plotit <- function(var, net, dyn, optcolor, randcolor, labelsize, df,
                   show.legend = "no", show.xlabel = "no", show.ylabel = "no") {
    labellookup <- switch(
        var,
        knn = expression(italic(k)[nn]),
        lcl = expression(italic(C)),#"C",
        pairs = expression(italic(K))#"K"
    )
    
    plotdf <- subset(df, dynamics == dyn & network == net)
    opt <- subset(plotdf, ns.type == "opt")[, var]
    rand <- subset(plotdf, ns.type == "rand")[, var]
    seedbreaks <- switch(var, 15, pairs = length(unique(plotdf$pairs)))
    
    hist.all <- hist(c(opt, rand), breaks = seedbreaks, plot = FALSE) # for the breaks
    hist.opt <- hist(opt, breaks = hist.all$breaks, plot = FALSE)
    hist.rand <- hist(rand, breaks = hist.all$breaks, plot = FALSE)

    ylim <- c(0, max(c(hist.opt$counts, hist.rand$counts))); print(ylim)

    plot(hist.opt, col = adjustcolor(optcolor, .5), ylim = ylim, main = "",
         xlab = switch(show.xlabel, no = "", yes = labellookup),
         ylab = switch(show.ylabel, no = "", yes = "Frequency"),
         cex.axis = labelsize, cex.lab = labelsize)
    plot(hist.rand, col = adjustcolor(randcolor, .5), add = TRUE)
    if(show.legend == "yes") {
        legend(
            "topright", cex = .75*labelsize, pch = 22, pt.bg = adjustcolor(c(optcolor, randcolor), .5), pt.cex = 2,
            legend = c("Optimized", "Random"), bty = "n"
        )
    }
}

ht <- wd <- 8
labelsize <- 2
optcolor <- "#3584e4"
randcolor <- "#33d17a"

if(save_plots) {
    pdf("../img/hist-knn.pdf", height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
par(mfcol = c(5, 2), mar = c(4, 4.5, 1.75, 0.5)) # keeping mar consistent makes the panels come out more evenly but hides the xlabels at the bottom. ylabels are also hidden. 
show.legend <- c("yes", rep("no", 9))
show.xlabel <- c(rep("no", 4), "yes", rep("no", 4), "yes")
show.ylabel <- c(rep("yes", 5), rep("no", 5))
for(i in seq_along(networks)) {
    plotit("knn", networks[i], "doublewell", optcolor, randcolor, labelsize, df,
           show.legend = show.legend[i], show.xlabel = show.xlabel[i], show.ylabel = show.ylabel[i])
    mtext(paste0("(", letters[i], ")"), cex = 0.75*labelsize, line = 0, adj = 0)
}
if(save_plots) dev.off()

if(save_plots) {
    pdf("../img/hist-lcl.pdf", height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
par(mfcol = c(5, 2), mar = c(4, 4.5, 1.9, 1.15))
show.legend <- c("yes", rep("no", 9))
show.xlabel <- c(rep("no", 4), "yes", rep("no", 4), "yes")
show.ylabel <- c(rep("yes", 5), rep("no", 5))
for(i in seq_along(networks)) {
    plotit("lcl", networks[i], "doublewell", optcolor, randcolor, labelsize, df,
           show.legend = show.legend[i], show.xlabel = show.xlabel[i], show.ylabel = show.ylabel[i])
    mtext(paste0("(", letters[i], ")"), cex = 0.75*labelsize, line = 0.1, adj = 0)
}
if(save_plots) dev.off()


morediscrete <- function(net, dyn, optcolor, randcolor, labelsize, df, nodesets,
                         show.legend = "no", show.xlabel = "no", show.ylabel = "no") {
    nbins <- length(nodesets[[paste(c(net, dyn), collapse = "_")]][[1]][[1]]$vs)

    plotdf <- subset(df, dynamics == dyn & network == net)
    opt <- subset(plotdf, ns.type == "opt")[, "pairs"]
    rand <- subset(plotdf, ns.type == "rand")[, "pairs"]

    table.opt <- tabulate(opt, nbins = nbins)
    table.rand <- tabulate(rand, nbins = nbins)

    ylim <- c(0, max(c(table.opt, table.rand)))

    barplot(
        matrix(c(table.opt, table.rand), byrow = TRUE, nrow = 2),
        col = adjustcolor(c(optcolor, randcolor), .5), beside = TRUE,
        ylim = ylim, cex.axis = labelsize,
        names.arg = seq(nbins), cex.names = labelsize, cex.lab = labelsize,
        ylab = switch(show.ylabel, no = "", yes = "Frequency"), #yaxt = "no",
        xlab = switch(show.xlabel, no = "", yes = expression(italic(K)))#"K")
    )
    if(show.legend == "yes") {
        legend(
            "topright", bty = "n", pch = 22, pt.bg = adjustcolor(c(optcolor, randcolor), .5), pt.cex = 2,
            legend = c("Optimized", "Random"), cex = .75*labelsize
        )
    }
}

if(save_plots) {
    pdf("../img/hist-comm.pdf", height = (ht/5)*3, width = wd)
} else {
    dev.new(height = (ht/5)*3, width = wd)
}
par(mfcol = c(3, 2), mar = c(4, 4.5, 1.9, 0.5))
show.legend <- c("yes", rep("no", 5))
show.xlabel <- rep(c(rep("no", 2), "yes"), 2)
show.ylabel <- c(rep("yes", 3), rep("no", 3))
commnets <- c(networks[1:5], "lfr")
for(i in seq_along(commnets)) {
    morediscrete(commnets[i], "doublewell", optcolor, randcolor, labelsize, df, nodesets,
                 show.legend = show.legend[i], show.xlabel = show.xlabel[i], show.ylabel = show.ylabel[i])
    mtext(paste0("(", letters[i], ")"), cex = 0.75*labelsize, line = 0.3, adj = 0)
}
if(save_plots) dev.off()






## Add new networks
## add k-core
## add betweenness and closeness, though these will be slow for large networks
networks <- c(
    "dolphin", "proximity", "celegans", "euroroad", "email",
    "drosophila", "reactome", "route_views", "spanish", "foldoc", "tree_of_life"
)
for(network in networks) {
    g <- readRDS(paste0("../data/", network, ".rds"))
    print(network)
    print(system.time(closeness(g)))
}

library(parallel)
ncores <- detectCores()-1
library(igraph)

set.seed(123) # because of the Louvain algorithm

save_plots <- FALSE # TRUE

useall <- "no" # "yes"    # If no, only opt and rand
useweights <- "no" # "yes"
weightsflag <- switch(useweights, no = FALSE, yes = TRUE)

get_from_ns <- function(row, ref, varname, collection) {
    network <- row[1]
    dynamic <- row[2]
    ns.type <- row[3]
    nslist <- paste(c(network, dynamic), collapse = "_")

    ns <- collection[[nslist]][[ns.type]]
    ref <- ref[[network]]

    var <- sapply(ns, function(x) mean(ref[x$vs], na.rm = TRUE))
    df <- data.frame(network = network, dynamics = dynamic, ns.type = ns.type, var = var,
                     row.names = NULL)
    colnames(df)[4] <- varname
    return(df)
}

max_fromsame <- function(ns, part) {
    mbr <- part[ns$vs]
    max(table(mbr))
}

get_pairs_scores <- function(row, collection, partitions) {
    network <- row[1]
    dynamics <- row[2]
    ns.type <- row[3]
    nslist <- paste(c(network, dynamics), collapse = "_")
    partition <- partitions[[network]]

    ns <- collection[[nslist]][[ns.type]]

    var <- sapply(ns, max_fromsame, partition)
    data.frame(network = network, dynamics = dynamics, ns.type = ns.type, pairs = var, row.names = NULL)
}


networks <- c( # try all networks
    "dolphin", "proximity", "celegans", "euroroad", "email", "er", "ba", "hk", "gkk", "lfr"
)
graphlist <- lapply(networks, function(network) readRDS(paste0("../data/", network, ".rds")))
dynamics <- c("doublewell", "SIS", "mutualistic", "genereg")
ns.types <- c("opt", "fixed", "rand", "constr", "quant", "knnconstr", "comm")
ns.types <- ns.types[switch(useall, no = c(1, 3), yes = 1:length(ns.types))]
conds <- expand.grid(networks, dynamics, ns.types, stringsAsFactors = FALSE)
colnames(conds) <- c("networks", "dynamics", "ns.type")
nslistnames <- apply(conds[, 1:2], 1, paste, collapse = "_")

fullstates <- lapply(networks, function(network) readRDS(paste0("../data/fullstate-", network, ".rds")))
names(fullstates) <- networks
nodesets <- lapply(nslistnames, function(cond) {
    readRDS(paste0("../data/ns-", cond, switch(useweights, no = ".rds", yes = "_w.rds")))
}) # _w if useweights
names(nodesets) <- nslistnames

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

#### NEW MODEL
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






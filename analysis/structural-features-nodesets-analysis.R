library(sfsmisc)

save_plots <- FALSE # TRUE

df <- read.csv("../data/nodeset-features.csv")

networks <- c(
    "dolphin", "proximity", "celegans", "euroroad", "email",
    "drosophila", "reactome", "route_views", "spanish", "foldoc",
    "tree_of_life", "word_assoc", "enron", "marker_cafe", "prosper",
    "er", "ba", "hk", "gkk", "lfr"
)
df$network <- factor(df$network, levels = networks)
df$dynamics <- factor(df$dynamics, levels = c("doublewell", "mutualistic", "SIS", "genereg"))
df$ns.type <- factor(df$ns.type, levels = c("rand", "opt", "fixed"))

ns <- sapply(networks, function(net) floor(log(igraph::vcount(readRDS(paste0("../data/", net, ".rds"))))))

                                        # Binary classification
model <- glm(
    ns.type ~ network + dynamics + k + cc + bc + knn + lcl + kcore + pairs + geods,
    family = binomial,
    data = subset(df, network != "er" & ns.type != "fixed")
)
summary(model)
print(1 - (model$deviance/model$null.deviance)) # 0.01640984

                                        # for tables?
exp(coefficients(model))
(exp(coefficients(model)) - 1)*100

checkit <- function(rng) {
    if(rng[1] > floor(rng[1]) & rng[1] < floor(rng[1]) + 0.25) rng[1] <- floor(rng[1])
    if(rng[2] < ceiling(rng[2]) & rng[2] > ceiling(rng[2]) - 0.25) rng[2] <- ceiling(rng[2])
    return(rng)
}

plotdistances <- function(net, dyn) {
    plotdf <- df[df$dynamics == dyn & df$network == net, ]
    opt <- plotdf$geods[plotdf$ns.type == "opt"]
    rand <- plotdf$geods[plotdf$ns.type == "rand"]

    hist.all <- hist(c(opt, rand), breaks = 15, plot = FALSE)
    hist.opt <- hist(opt, breaks = hist.all$breaks, plot = FALSE)
    hist.rand <- hist(rand, breaks = hist.all$breaks, plot = FALSE)
    ylim <- range(c(hist.opt$counts, hist.rand$counts))

    xlim <- checkit(range(hist.all$breaks))
    ## if(xlim[1] > floor(xlim[1]) & xlim[1] < floor(xlim[1]) + 0.25) xlim[1] <- floor(xlim[1])
    ## if(xlim[2] < ceiling(xlim[2]) & xlim[2] > ceiling(xlim[2]) - 0.25) xlim[2] <- ceiling(xlim[2])

    plot(hist.rand, col = adjustcolor(randcolor, .5), ylim = ylim, xlim = xlim, main = "",
         xlab = "", ylab = "", axes = FALSE)
    eaxis(1, cex.axis = 1.5)
    eaxis(2, cex.axis = 1.5)
    plot(hist.opt, col = adjustcolor(optcolor, .5), add = TRUE)
}

wd <- 12
ht <- 10
optcolor <- "#3584e4"
randcolor <- "#33d17a"
if(save_plots) {
    pdf("../img/hist-dist-v2.pdf", height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
par(mfrow = c(4, 5), mar = c(2, 2.5, 0.5, 0.5), pty = "s")
for(i in seq_along(networks)) {
    plotdistances(networks[i], "doublewell")
    mtext(paste0("(", letters[i], ")"), cex = 1.5, line = -2, adj = 0.02)
    if(i == 1) {
        legend("left", legend = c("Optimized", "Random"), bty = "n",
               pch = 22, pt.bg = adjustcolor(c(optcolor, randcolor), .5), pt.cex = 2)
    }
}
if(save_plots) dev.off()

plotpairs <- function(net, dyn, labelsize = 1.5) {
    plotdf <- df[df$dynamics == dyn & df$network == net, ]
    opt <- plotdf$pairs[plotdf$ns.type == "opt"]
    rand <- plotdf$pairs[plotdf$ns.type == "rand"]

    nbins <- ns[net]
    table.opt <- tabulate(opt, nbins = nbins)
    table.rand <- tabulate(rand, nbins = nbins)

    ylim <- c(0, max(c(table.opt, table.rand)))
    
    barplot(
        matrix(c(table.opt, table.rand), byrow = TRUE, nrow = 2),
        col = adjustcolor(c(optcolor, randcolor), .5), beside = TRUE,
        ylim = ylim, cex.axis = labelsize, ylab = "", xlab = "", 
        names.arg = seq(nbins), cex.names = labelsize, cex.lab = labelsize, axes = FALSE
    )
    eaxis(2, cex.axis = labelsize)
}

wd <- 12
ht <- 10
optcolor <- "#3584e4"
randcolor <- "#33d17a"
if(save_plots) {
    pdf("../img/hist-comm-v2.pdf", height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
par(mfrow = c(4, 4), mar = c(2, 2.5, 0.5, 0.5), pty = "s")
for(i in seq_along(networks[-which(networks %in% c("er", "ba", "hk", "gkk"))])) {
    plotpairs(networks[i], "doublewell")
    mtext(paste0("(", letters[i], ")"), cex = 1.5, line = -2, adj = 0.02)
    if(i == 1) {
        legend("topright", legend = c("Optimized", "Random"), bty = "n",
               pch = 22, pt.bg = adjustcolor(c(optcolor, randcolor), .5), pt.cex = 2)
    }
}
if(save_plots) dev.off()


plotks <- function(net, dyn) {
    ## net <- "route_views"
    plotdf <- df[df$dynamics == dyn & df$network == net, ]
    opt <- plotdf$k[plotdf$ns.type == "opt"]
    ## fixed <- plotdf$k[plotdf$ns.type == "fixed"]
    rand <- plotdf$k[plotdf$ns.type == "rand"]

    hist.all <- hist(c(opt, rand), breaks = 15, plot = FALSE)
    hist.opt <- hist(opt, breaks = hist.all$breaks, plot = FALSE)
    hist.rand <- hist(rand, breaks = hist.all$breaks, plot = FALSE)

    ylim <- range(c(hist.opt$counts, hist.rand$counts))
    xlim <- range(hist.all$breaks)
    if(net %in% networks[which(letters %in% c("f", "h", "l"))]) xlim[1] <- 0

    plot(hist.rand, col = adjustcolor(randcolor, .5), ylim = ylim, xlim = xlim, main = "",
         xlab = "", ylab = "", axes = FALSE)
    eaxis(1, cex.axis = 1.5)
    eaxis(2, cex.axis = 1.5)
    plot(hist.opt, col = adjustcolor(optcolor, .5), add = TRUE)
}

wd <- 12
ht <- 10
optcolor <- "#3584e4"
## fixedcolor <- "#ff7800"
randcolor <- "#33d17a"
if(save_plots) {
    pdf("../img/hist-ks-v2.pdf", height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
par(mfrow = c(4, 5), mar = c(2, 2.5, 0.5, 0.5), pty = "s")
for(i in seq_along(networks)) {
    plotks(networks[i], "doublewell")
    mtext(paste0("(", letters[i], ")"), cex = 1.5, line = -2, adj = 0.02)
    if(i == 1) {
        legend("left", legend = c("Optimized", "Random"), bty = "n",
               pch = 22, pt.bg = adjustcolor(c(optcolor, randcolor), .5), pt.cex = 2)
    }
}
if(save_plots) dev.off()


scattereps <- function(net, dyn) {
    plotdf <- df[df$dynamics == "doublewell" & df$network == net, ]
    xlim <- range(plotdf$k)
    ylim <- range(plotdf$error)
    
    plot(NULL, xlim = xlim, ylim = ylim, xlab = "", ylab = "", axes = FALSE, log = "y")
    points(error ~ k, data = plotdf, subset = ns.type == "rand", pch = 0, col = randcolor)
    points(error ~ k, data = plotdf, subset = ns.type == "fixed", pch = 2, col = fixedcolor)
    points(error ~ k, data = plotdf, subset = ns.type == "opt", pch = 1, col = optcolor)
    box()
    eaxis(1, cex.axis = 1.5)
    eaxis(2, cex.axis = 1.5, n.axp = 1)
}


wd <- 12
ht <- 10
optcolor <- "#3584e4"
fixedcolor <- "#ff7800"
randcolor <- "#33d17a"
if(save_plots) {
    pdf("../img/k-scatter.pdf", height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
par(mfrow = c(4, 5), mar = c(2, 3.5, 0.5, 0.5), pty = "s")
for(i in seq_along(networks)) {
    scattereps(networks[i], "doublewell")
    mtext(paste0("(", letters[i], ")"), cex = 1.5, line = -2, adj = 0.02)
    if(i == 1) {
        legend("topright", legend = c("Optimized", "Degree-preserving", "Random"), bty = "n",
               pch = c(1, 2, 0), col = c(optcolor, fixedcolor, randcolor), pt.cex = 1.5, pt.lwd = 2)
    }
}
if(save_plots) dev.off()

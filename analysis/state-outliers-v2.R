library(igraph)
library(optNS)
library(sfsmisc)

save_plots <- FALSE # TRUE

get_zs <- function(vs, fs, dispersion = c("mad", "sd"), center = c("mean", "median")) { 
    dispersion <- get(dispersion)
    center <- get(center)

    m <- center(fs)
    s <- dispersion(fs)
    if(s == 0) {
        return(rep(NA, length(vs)))
    } else {
        (fs[vs] - m)/s
    }
}

nets <- c(
    "dolphin", "proximity", "celegans", "euroroad", "email",
    "drosophila", "reactome", "route_views", "spanish", "foldoc",
    "tree_of_life", "word_assoc", "enron", "marker_cafe", "prosper",
    "er", "ba", "hk", "gkk", "lfr"
)
## nets <- nets[-which(nets == "er")]
dynamics <- c(
    "doublewell", "mutualistic", "SIS", "genereg"
)

net <- "proximity" # dolphin
dispersion <- "mad" # sd
center <- "median" # mean

fss <- readRDS(paste0("../data/fullstate-", net, ".rds"))

## only show z values, not x values, so convert up front.
analyze <- function(net, dyn, dispersion, center) {
    fs <- fss[[dyn]]
    nss <- readRDS(paste0("../data/ns-", net, "_", dyn, ".rds"))$opt

    fsz <- t(apply(fs, 1, function(xi) get_zs(seq_len(length(xi)), xi, dispersion, center)))

    Ds <- switch(dyn, seq(0, 1, length.out = 100), mutualistic = seq(0, 3, length.out = 100))

    matplot(Ds, fsz, type = "l", lty = 1, lwd = 0.5, col = adjustcolor(8, 0.4),
            xlab = "D", ylab = "z", font.lab = 3, cex.lab = 2, axes = FALSE)
    for(ns in nss) {
        matlines(Ds, fsz[, ns$vs], lty = 1, lwd = 0.5, col = "#3584e4")
    }
    box()
    eaxis(1, cex.axis = 1.5)
    eaxis(2, cex.axis = 1.5)

    par(new = TRUE)
    disps <- apply(fs, 1, get(dispersion))
    dispcolor <- 3
    plot(Ds, disps, type = "l", lty = 1, lwd = 2, col = dispcolor, axes = FALSE, xlab = "", ylab = "")
    eaxis(4, cex.axis = 1.5, col.axis = dispcolor, col = dispcolor, small.args = list(col = dispcolor))
}

ht <- 10
wd <- 10
if(save_plots) {
    pdf(paste0("../img/state-outliers-", dispersion, "-", net, "-zonly.pdf"), height = ht, width = wd)
} else {
    dev.new(width = wd, height = ht)
}
par(mfrow = c(2, 2), mar = c(5, 5, 1, 4), pty = "s")
for(i in seq_along(dynamics)) {
    analyze(net, dynamics[i], dispersion, center)
    if(i == 1) {
        legend("bottomright", bty = "n",
               legend = c(
                   "Node state",
                   "Sentinel node state",
                   "Measure of dispersion"
               ),
               col = c(8, "#3584e4", 3), lty = 1, lwd = 1.5)
    }
    mtext(dynamics[i], line = -2, cex = 1.5)
}
if(save_plots) dev.off()




## Two D values above the transition. Say, the first one above the transition and the final. I should have a way to determine the former?
## At each D, normalized ([0, 1]?) histogram for the sentinel node z-values and all node z-values.
## Count fraction of sentinel nodes attaining |z| > θ, where θ = 2 or some other reasonable value. 1.96 is 2σ.
## Do all of that for four dynamics, two networks. What we have above.
##
## Let's start with the final for now.
conds <- expand.grid(nets, dynamics, stringsAsFactors = FALSE)
colnames(conds) <- c("network", "dynamics")

fss <- lapply(nets, function(net) readRDS(paste0("../data/fullstate-", net, ".rds")))
names(fss) <- nets
nss <- apply(conds, 1, function(row) {
    readRDS(paste0("../data/ns-", row["network"], "_", row["dynamics"], ".rds"))$opt
}, simplify = FALSE)
names(nss) <- apply(conds, 1, paste, collapse = "_")

zanalysis <- function(net, dyn, ell = 100) {
    fs <- fss[[net]][[dyn]]
    ns <- nss[[paste(net, dyn, sep = "_")]]
    xi <- fs[ell, ]
    zi <- get_zs(seq_len(length(xi)), xi, "mad", "median")
    theta_zi <- quantile(abs(zi), 0.95, names = FALSE)

    vs <- as.numeric(get_vs(ns))
    vsz <- get_zs(vs, xi, "mad", "median")
    theta_vs <- quantile(abs(vsz), 0.95, names = FALSE)
                                        # now the z value representing 95% of |z|
    return(data.frame(all = theta_zi, sentinels = theta_vs, net = net, dyn = dyn, N = length(xi))) 
}

res <- data.frame(
    do.call(
        rbind,
        apply(conds, 1, function(row) zanalysis(row["network"], row["dynamics"]))
    )
)

                                        # positive t means first arg is larger than second arg
with(list(df = subset(res, dyn == "doublewell")), t.test(df$all, df$sentinels, paired = TRUE))
                                        # t = 2.0692, df = 19, p-value = 0.05241
with(list(df = subset(res, dyn == "mutualistic")), t.test(df$all, df$sentinels, paired = TRUE))
                                        # t = 2.056, df = 19, p-value = 0.05379
with(list(df = subset(res, dyn == "SIS")), t.test(df$all, df$sentinels, paired = TRUE))
                                        # t = 2.1041, df = 19, p-value = 0.04891
with(list(df = subset(res, dyn == "genereg")), t.test(df$all, df$sentinels, paired = TRUE))
                                        # t = -1.0627, df = 19, p-value = 0.3012

ht <- wd <- 10
dev.new(height = ht, width = wd)
par(mfrow = c(2, 2), mar = c(5, 5, 1, 1), pty = "s")
for(dynamic in dynamics) {
    df <- subset(res, dyn == dynamic)
    print(head(df))
    matplot(
        df$N, df[, c("all", "sentinels")], type = "p", pch = c(0, 1), col = 1:2, cex = 2, lwd = 2, log = "x",
        xlab = "", ylab = "", axes = FALSE
    )
    box()
    eaxis(1, n.axp = 1, cex.axis = 1.5)
    eaxis(2, cex.axis = 1.5)
    title(xlab = "N", cex.lab = 2, font.lab = 3)
    title(ylab = "95th percentile of |z|", cex.lab = 2)#, line = 5)
    if(dynamic == "doublewell") {
        legend(
            "topleft", bty = "n", legend = c("All", "Sentinel"),
            col = 1:2, pch = 0:1, pt.cex = 1.5, pt.lwd = 1.5
        )
    }
    title(main = dynamic)
}



test <- subset(res, dyn == "doublewell")


res <- as.data.frame(t(sapply(nets, function(net) {
    ## net <- "ba"
    dyn <- "doublewell" # genereg
    ## theta <- qnorm(0.975)
    ell <- 100 # the index on Ds
    fs <- readRDS(paste0("../data/fullstate-", net, ".rds"))[[dyn]]
    nss <- readRDS(paste0("../data/ns-", net, "_", dyn, ".rds"))$opt

    xi <- fs[ell, ]
    zi <- get_zs(seq_len(length(xi)), xi, "mad", "median")
    theta_zi <- quantile(abs(zi), 0.95, names = FALSE)

    vs <- as.numeric(get_vs(nss))
    vsz <- get_zs(vs, xi, "mad", "median")
    theta_vs <- quantile(abs(vsz), 0.95, names = FALSE)
    ## return(c(all = sum(abs(zi) >= theta)/length(zi), sentinels = sum(abs(vsz) >= theta)/length(vsz)))
    return(c(all = theta_zi, sentinels = theta_vs)) # now the z value representing 95% of |z|
})))
t.test(res$all, res$sentinels, paired = TRUE)
res$Ns <- sapply(nets, function(net) vcount(readRDS(paste0("../data/", net, ".rds"))))



hists <- list(all = hist(zi, breaks = 20, plot = FALSE))
hists$sentinels <- hist(vsz, breaks = hists$all$breaks, plot = FALSE)

hists$all$counts <- hists$all$counts/sum(hists$all$counts)
hists$sentinels$counts <- hists$sentinels$counts/sum(hists$sentinels$counts)

ylim <- range(c(hists$all$counts, hists$sentinels$counts))

dev.new()
par(mar = c(5, 6, 1, 1), pty = "s")
plot(hists$all, col = adjustcolor(1, .5), ylim = ylim, axes = FALSE, xlab = "", ylab = "", main = "")
plot(hists$sentinels, col = adjustcolor(2, .5), add = TRUE, ylim = ylim)
box()
eaxis(1, cex.axis = 1.5)
eaxis(2, cex.axis = 1.5)
title(xlab = "z", cex.lab = 2, font.lab = 3)
title(ylab = "Normalized frequency", cex.lab = 2, line = 4)
title(main = dyn)


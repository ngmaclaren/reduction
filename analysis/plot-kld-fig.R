library(sfsmisc)

save_plots <- TRUE # FALSE

CI <- function(x, z) mean(x) + z*(sd(x)/sqrt(length(x)))

load("kld-fig.RData")
ht <- 5
wd <- 10
labelsize <- 1.75

if(save_plots) {
    pdf("../img/kld-fig.pdf", height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}

par(mar = c(4.5, 5.5, 1.5, 0.5), mfrow = c(1, 2))
xlim <- c(1, 12)
optcol <- "#3584e4"
optpch <- 16 
randcol <- "#33d17a"
randpch <- 15

with(dolphin.env, {
    ylim <- range(c(as.numeric(KLs.rands), KLs.opts))
    plot(NULL, xlim = xlim, ylim = ylim, xlab = "", ylab = "", xaxt = "n", yaxt = "n", log = "y")
    eaxis(1, cex.axis = labelsize)
    eaxis(2, n.axp = 1, cex.axis = labelsize)
    title(xlab = expression(italic(n)),
          cex.lab = labelsize)
    title(ylab = "Kullback-Leibler divergence", line = 4, cex.lab = labelsize)
    mtext("(a)", line = 0.2, adj = -0.1, cex = labelsize)
    points(seq(12), KLs.opts, col = optcol, pch = optpch, cex = 3)
    points(seq(12), colMeans(KLs.rands), col = randcol, pch = randpch, cex = 3)
    segments(
        x0 = seq(12),
        y0 = apply(KLs.rands, 2, CI, -1.96),
        y1 = apply(KLs.rands, 2, CI, 1.96),
        col = randcol, lwd = 2
    )
    legend(
        "topright", bty = "n", col = c(optcol, randcol), pch = c(optpch, randpch), pt.cex = 2, cex = 0.75*labelsize,
        legend = c("Optimized", "Random")
    )
})

with(ba.env, {
    ylim <- range(c(as.numeric(KLs.rands), KLs.opts))
    plot(NULL, xlim = xlim, ylim = ylim, xlab = "", ylab = "", xaxt = "n", yaxt = "n", log = "y")
    eaxis(1, cex.axis = labelsize)
    eaxis(2, n.axp = 1, cex.axis = labelsize)
    title(xlab = expression(italic(n)),#"n",
          cex.lab = labelsize)
    title(ylab = "Kullback-Leibler divergence", line = 4, cex.lab = labelsize)
    mtext("(b)", line = 0.2, adj = -0.1, cex = labelsize)
    points(seq(12), KLs.opts, col = optcol, pch = optpch, cex = 3)
    points(seq(12), colMeans(KLs.rands), col = randcol, pch = randpch, cex = 3)
    segments(
        x0 = seq(12),
        y0 = apply(KLs.rands, 2, CI, -1.96),
        y1 = apply(KLs.rands, 2, CI, 1.96),
        col = randcol, lwd = 2
    )
})
       
if(save_plots) dev.off()

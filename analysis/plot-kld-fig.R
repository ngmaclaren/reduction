library(sfsmisc)

save_plots <- TRUE # FALSE

load("kld-fig.RData")
ht <- 7
wd <- 14
labelsize <- 1.75
palette("Tableau 10")

if(save_plots) {
    pdf("../img/kld-fig.pdf", height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}

par(mar = c(5, 5.5, 1, 0.5), mfrow = c(1, 2))
xlim <- c(1, 12)

ylim <- range(c(as.numeric(KLs.rands), KLs.opts))
plot(NULL, xlim = xlim, ylim = ylim, xlab = "", ylab = "", xaxt = "n", yaxt = "n", log = "y")
eaxis(1, cex.axis = labelsize)
eaxis(2, n.axp = 1, cex.axis = labelsize)
title(xlab = "n", cex.lab = labelsize)
title(ylab = "Kullback-Leibler divergence", line = 4, cex.lab = labelsize) 
points(seq(12), KLs.opts, col = 1, pch = 16, cex = 3)
points(seq(12), colMeans(KLs.rands), col = 2, pch = 15, cex = 3)
segments(
    x0 = seq(12),
    y0 = apply(KLs.rands, 2, quantile, probs = 0.025),
    y1 = apply(KLs.rands, 2, quantile, probs = 0.975),
    col = 2, lwd = 2
)

ylim <- range(c(as.numeric(errors.rands), as.numeric(errors.opts)))
plot(NULL, xlim = xlim, ylim = ylim, xlab = "", ylab = "", xaxt = "n", yaxt = "n", log = "y")
eaxis(1, cex.axis = labelsize)
eaxis(2, n.axp = 1, cex.axis = labelsize)
title(xlab = "n", cex.lab = labelsize)
title(ylab = "Approximation error", line = 4, cex.lab = labelsize)
points(seq(12), colMeans(errors.opts), col = 1, pch = 16, cex = 3)
points(seq(12), colMeans(errors.rands), col = 2, pch = 15, cex = 3)
segments(
    x0 = seq(12),
    y0 = apply(errors.rands, 2, quantile, probs = 0.025),
    y1 = apply(errors.rands, 2, quantile, probs = 0.975),
    col = 2, lwd = 2
)

if(save_plots) dev.off()

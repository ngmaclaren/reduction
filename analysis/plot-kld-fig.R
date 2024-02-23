load("kld-fig.RData")

ht <- 7
wd <- 14
labelsize <- 1.75
palette("Tableau 10")

dev.new(height = ht, width = wd)
par(mar = c(5, 5, 1, 1), mfrow = c(1, 2))
xlim <- c(1, 12)
ylim <- range(c(as.numeric(KLs.rands), KLs.opts))
plot(NULL, xlim = xlim, ylim = ylim, xlab = "n", ylab = "Kullback-Leibler divergence", log = "y",
     cex.axis = labelsize, cex.lab = labelsize)
points(seq(12), KLs.opts, col = 1, pch = 16, cex = 3)
points(seq(12), colMeans(KLs.rands), col = 2, pch = 16, cex = 3)
segments(
    x0 = seq(12),
    y0 = apply(KLs.rands, 2, quantile, probs = 0.025),
    y1 = apply(KLs.rands, 2, quantile, probs = 0.975),
    col = 2, lwd = 2
)
ylim <- range(c(as.numeric(errors.rands), as.numeric(errors.opts)))
plot(NULL, xlim = xlim, ylim = ylim, xlab = "n", ylab = "Approximation error", log = "y",
     cex.axis = labelsize, cex.lab = labelsize)
points(seq(12), colMeans(errors.opts), col = 1, pch = 16, cex = 3)
points(seq(12), colMeans(errors.rands), col = 2, pch = 16, cex = 3)
segments(
    x0 = seq(12),
    y0 = apply(errors.rands, 2, quantile, probs = 0.025),
    y1 = apply(errors.rands, 2, quantile, probs = 0.975),
    col = 2, lwd = 2
)

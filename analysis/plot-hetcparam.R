library(sfsmisc)

network <- "proximity"
save_plots <- TRUE # FALSE

load(paste0("heterogeneous-cparam-", network, ".RData"))

optcolor <- "#3584e4"
randcolor <- "#33d17a"
grays <- gray(probs)
adjust <- seq(-0.5, 0.5, length.out = length(probs) + 2)
adjust <- adjust[c(-1, -length(adjust))]
ptcex <- 0.9
amt <- 0.015

wd <- 8
if(save_plots) {
    pdf(paste0("../img/heterogeneous-cparam-", network, ".pdf"), width = wd)
} else {
    dev.new(width = wd)
}
par(mar = c(5, 16, 1, 1))#, pty = "s")
plot(NULL, ylim = rev(c(0.85, 6.4)), xlim = xlim, axes = FALSE, xlab = "", ylab = "", log = "x")
points(rand_error, jitter(rep(1, 100), amount = 0.1), col = randcolor, pch = 0)
points(opt_error, jitter(rep(1, 100), amount = 0.1), col = optcolor, pch = 1)
for(i in seq_along(Ylist)) {
    adj <- adjust[i]
    points(
        het_errors[[i]]$random, jitter(rep(6 + adj, 100), amount = amt), col = grays[i], pch = 1, cex = ptcex
    )
    points(# Error in xy.coords(x, y) : 'x' and 'y' lengths differ
        het_errors[[i]]$degsum.high, jitter(rep(5 + adj, 100), amount = amt), col = grays[i], pch = 1, cex = ptcex
    )
    points(
        het_errors[[i]]$degsum.low, jitter(rep(4 + adj, 100), amount = amt), col = grays[i], pch = 1, cex = ptcex
    )
    points(
        het_errors[[i]]$edgebet.high, jitter(rep(3 + adj, 100), amount = amt), col = grays[i], pch = 1, cex = ptcex
    )
    points(
        het_errors[[i]]$edgebet.low, jitter(rep(2 + adj, 100), amount = amt), col = grays[i], pch = 1, cex = ptcex
    )
}
box()
eaxis(1, cex.axis = 1.5, n.axp = 1)
title(xlab = "Approximation error", cex.lab = 2, line = 3.5)
axis(
    2, at = 1:6, tick = FALSE, las = 2, cex.axis = 1.5,
    labels = c(
        "Original",
        "Lowest edge betweenness", #"Edge betweenness: low"
        "Highest edge betweenness", #"Edge betweenness: high",
        "Lowest degree sum", #"Degree sum: low",
        "Highest degree sum", #"Degree sum: high",
        "Random"
    )
)
if(save_plots) dev.off()

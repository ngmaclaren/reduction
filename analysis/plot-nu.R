library(sfsmisc)
library(optNS)

network <- "proximity"
load(paste0("nu-", network, ".RData"))
save_plots <- FALSE # TRUE

rand_error <- get_error(nodesets$rand)
opt_error <- get_error(nodesets$opt)

randcolor <- "#33d17a"
optcolor <- "#3584e4"

## EDGES

ylim <- range(c(unlist(uncertain_errorE), rand_error, opt_error))

ht <- 5
wd <- 8
if(save_plots) {
    pdf(paste0("../img/edge-deletion-", network, ".pdf"), width = wd, height = ht)
} else {
    dev.new(width = 8, height = ht)
}
lyt <- matrix(c(
    1, 1, 1, 2,
    1, 1, 1, 2), byrow = TRUE, nrow = 2)
layout(lyt)
par(mar = c(5, 6, 1, 0))
matplot(edges_toremove/M, uncertain_errorE, log = "y",
     pch = 1, col = 1, cex = 2, lwd = 2, ylim = ylim,
     axes = FALSE, xlab = "", ylab = "")
box()
eaxis(1, cex.axis = 1.5)
eaxis(2, cex.axis = 1.5, n.axp = 1)
title(xlab = "Proportion of edges removed", cex.lab = 2)
title(ylab = "Approximation error", cex.lab = 2, line = 3.75)
legend("topleft", bty = "n", pch = c(1, 0, 1), col = c(optcolor, randcolor, 1), cex = 1.5,
       pt.lwd = 2, pt.cex = 2,
       legend = c("Optimized, full information", "Random", "Optimized, partial information"))
par(mar = c(5, 0, 1, 1))
plot(
    jitter(rep(1, 100), amount = 0.1), rand_error, log = "y", ylim = ylim, xlim = c(0.85, 1.15),
    cex = 2, lwd = 2, pch = 0, col = randcolor, axes = FALSE, xlab = "", ylab = ""
)
points(jitter(rep(1, 100), amount = 0.1), opt_error, cex = 2, pch = 1, col = optcolor, lwd = 2)
if(save_plots) dev.off()

## NODES

ylim <- range(c(as.numeric(uncertain_errorN), rand_error, opt_error))
ht <- 5
wd <- 8
if(save_plots) {
    pdf(paste0("../img/node-deletion-", network, ".pdf"), width = wd, height = ht)
} else {
    dev.new(width = 8, height = ht)
}
lyt <- matrix(c(
    1, 1, 1, 2,
    1, 1, 1, 2), byrow = TRUE, nrow = 2)
layout(lyt)
par(mar = c(5, 6, 1, 0))
matplot((N - nodes_tokeep)/N, uncertain_errorN, log = "y",
     pch = 1, col = 1, cex = 2, lwd = 2, ylim = ylim,
     axes = FALSE, xlab = "", ylab = "")
box()
eaxis(1, cex.axis = 1.5)
eaxis(4, cex.axis = 1.5, n.axp = 1)
title(xlab = "Proportion of nodes removed", cex.lab = 2)
title(ylab = "Approximation error", cex.lab = 2, line = 3.75)
legend("topright", bty = "n", pch = c(1, 0, 1), col = c(optcolor, randcolor, 1), cex = 1.5,
       pt.lwd = 2, pt.cex = 2,
       legend = c("Optimized, full info", "Random", "Optimized, partial info"))
par(mar = c(5, 0, 1, 1))
plot(
    jitter(rep(1, 100), amount = 0.1), rand_error, log = "y", xlim = c(0.75, 1.25), ylim = ylim,
    cex = 2, lwd = 2, pch = 0,
    col = randcolor, axes = FALSE, xlab = "", ylab = ""
)
points(jitter(rep(1, 100), amount = 0.1), opt_error, cex = 2, pch = 1, col = optcolor, lwd = 2)
if(save_plots) dev.off()

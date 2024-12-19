dynamics <- "mutualistic" # genereg
infile <- paste0("SI-comps-", dynamics, ".RData")

library(optNS)
load(infile)
save_plots <- TRUE # FALSE
imgfile <- paste0("../img/SI-comps-", dynamics, ".pdf")
legendloc <- switch(dynamics, mutualistic = "topleft", genereg = "topleft")
labelsize <- 2.5
legendlabelsize <- 0.6*labelsize

                                        # opt
get_error(solns)

                                        # GBB/Cheby
calc_obj(GBB, y)
calc_obj(Cheby2, y)
calc_obj(Cheby3, y)

calc_obj(GBB, GBB.obs)
calc_obj(Cheby2, GBB.obs)
calc_obj(Cheby3, GBB.obs)

                                        # DART
apply(DARTs, 2, calc_obj, y = y)

sapply(1:3, function(i) calc_obj(DARTs[, i], DART.obs[, i]))

## for plotting, it's just the six panels of what I call bifurcation diagrams (but I'm not sure that's what they're actually called.

Ds <- params$Ds

ht <- 13*2/3
wd <- 13
if(save_plots) {
    pdf(imgfile, height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
par(mfrow = c(2, 3), mar = c(5, 5, 1, 1), pty = "s")
## a: SA, n = 1
plot_ns(Ds, Y, font.lab = 3)
lines(Ds, y, lwd = 6, col = "#000000")
lines(Ds, Y[, solns[[1]]$vs], lwd = 6, col = "#9141ac")
lines(Ds, Y[, solns[[1]]$vs], lwd = 4, col = "#ff7800")
mtext("(a)", line = -3, at = max(Ds), adj = 1, cex = 2.25)
legend(
    legendloc, bty = "n", lwd = c(2, 6, 6, 4), lty = 1, cex = legendlabelsize,
    col = c("#babdb6", "#000000", "#9141ac", "#ff7800"),
    legend = c("Node state", "Average state", "Sentinel node approximation", "Sentinel nodes, n = 1")
)
## b: SA, n = 2
plot_ns(Ds, Y, font.lab = 3)
lines(Ds, y, lwd = 6, col = "#000000")
matlines(Ds, Y[, solns[[2]]$vs], lty = 1, lwd = 4, col = "#f6d32d")
lines(Ds, rowMeans(Y[, solns[[2]]$vs]), lwd = 6, col = "#9141ac")
mtext("(b)", line = -3, at = max(Ds), adj = 1, cex = 2.25)
legend(
    legendloc, bty = "n", lwd = 4, lty = 1, cex = legendlabelsize,
    col = "#f6d32d", legend = "Sentinel nodes, n = 2"
)
## c: SA, n = 3
plot_ns(Ds, Y, font.lab = 3)
lines(Ds, y, lwd = 6, col = "#000000")
matlines(Ds, Y[, solns[[3]]$vs], lty = 1, lwd = 4, col = "#33d17a")
lines(Ds, rowMeans(Y[, solns[[3]]$vs]), lwd = 6, col = "#9141ac")
mtext("(c)", line = -3, at = max(Ds), adj = 1, cex = 2.25)
legend(
    legendloc, bty = "n", lwd = 4, lty = 1, cex = legendlabelsize,
    col = "#33d17a", legend = "Sentinel nodes, n = 3"
)
## d: SA, n = 4
plot_ns(Ds, Y, font.lab = 3)
lines(Ds, y, lwd = 6, col = "#000000")
matlines(Ds, Y[, solns[[4]]$vs], lty = 1, lwd = 4, col = "#3584e4")
lines(Ds, rowMeans(Y[, solns[[4]]$vs]), lwd = 6, col = "#9141ac")
mtext("(d)", line = -3, at = max(Ds), adj = 1, cex = 2.25)
legend(
    legendloc, bty = "n", lwd = 4, lty = 1, cex = legendlabelsize,
    col = "#3584e4", legend = "Sentinel nodes, n = 4"
)
## e: GBB, including original and Chebyshev order 2 and 3
plot_ns(Ds, Y, font.lab = 3)
lines(Ds, y, lwd = 6, col = "#000000")
lines(Ds, GBB.obs, lwd = 6, col = "#986a44")
lines(Ds, GBB, lwd = 6, col = "#e01b24")
lines(Ds, Cheby2, lwd = 4, col = "#e01b24") # , lty = 2
lines(Ds, Cheby3, lwd = 2, col = "#e01b24") # , lty = 3
mtext("(e)", line = -3, at = max(Ds), adj = 1, cex = 2.25)
legend(
    legendloc, bty = "n", cex = legendlabelsize,
    lty = 1, #c(1, 1, 2, 3),
    lwd = c(6, 6, 4, 2),
    col = c("#986a44", rep("#e01b24", 3)),
    legend = c(
        "GBB observable",
        "GBB approximation",
        "Chebyshev, order 2",
        "Chebyshev, order 3"
    )
)
## f: DART, including original and DART dim 2 and 3
plot_ns(Ds, Y, font.lab = 3)
lines(Ds, y, lwd = 6, col = "#000000")
lines(Ds, DART.obs[, 1], lwd = 6, col = "#e5a50a")
##lines(Ds, DART.obs[, 2], lwd = 6, col = "#f6d32d")
##lines(Ds, DART.obs[, 3], lwd = 6, col = "#e5a50a")
lines(Ds, DARTs[, 1], lwd = 6, col = "#ff007f")
lines(Ds, DARTs[, 2], lwd = 4, col = "#ff007f")
lines(Ds, DARTs[, 3], lwd = 2, col = "#ff007f")
mtext("(f)", line = -3, at = max(Ds), adj = 1, cex = 2.25)
legend(
    legendloc, bty = "n", cex = legendlabelsize,
    lty = 1,#c(1, 1, 2, 3),
    lwd = c(6, 6, 4, 2),
    col = c("#e5a50a", rep("#ff007f", 3)),
    legend = c(
        "DART observable",#, 1-D",
        ##"DART observable, 2-D",
        ##"DART observable, 3-D",
        "DART approx., 1-dim",
        "DART approx., 2-dim",
        "DART approx., 3-dim"
    )
)
if(save_plots) dev.off()

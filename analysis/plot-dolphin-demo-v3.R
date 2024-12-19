library(optNS)
load("dolphin-demo.RData")

save_plots <- TRUE # FALSE
imgfile <- "../img/dolphin-demo-v3.pdf"
labelsize <- 2.5
legendlabelsize <- 0.6*labelsize
legendloc <- "bottomright"

plot_ns <- function(Ds, Y, ...) {
    matplot(
        Ds, Y, type = "l", lty = 1, lwd = 0.5, col = "#babdb6",
        xlab = "D", ylab = "x", cex.lab = labelsize, cex.axis = labelsize, ...
    )
}


                                        # opt
get_error(solns)

                                        # GBB
calc_obj(GBB, y) # against xbar
calc_obj(GBB, GBB.obs) # against GBB observable

                                        # DART
calc_obj(DART, y) # against xbar
calc_obj(DART, DART.obs) # against DART observable

Ds <- sdn::.doublewell$Ds

ht <- 13*2/3
wd <- 13
if(save_plots) {
    pdf(imgfile, height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
par(mfrow = c(2, 3), mar = c(5, 5, 1, 1), pty = "s") # mar = c(5, 5, 1, 1) , mai = c(0.6, 0.6, 0.1, 0.1)
## a: SA, n = 1
plot_ns(Ds, Y, font.lab = 3)
lines(Ds, y, lwd = 6, col = "#000000")
lines(Ds, Y[, solns[[1]]$vs], lwd = 6, col = "#9141ac")
lines(Ds, Y[, solns[[1]]$vs], lwd = 4, col = "#ff7800")
mtext("(a)", line = -3, at = min(Ds), adj = 0, cex = 2.25)
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
mtext("(b)", line = -3, at = min(Ds), adj = 0, cex = 2.25)
legend(
    legendloc, bty = "n", lwd = 4, lty = 1, cex = legendlabelsize,
    col = "#f6d32d", legend = "Sentinel nodes, n = 2"
)
## c: SA, n = 3
plot_ns(Ds, Y, font.lab = 3)
lines(Ds, y, lwd = 6, col = "#000000")
matlines(Ds, Y[, solns[[3]]$vs], lty = 1, lwd = 4, col = "#33d17a")
lines(Ds, rowMeans(Y[, solns[[3]]$vs]), lwd = 6, col = "#9141ac")
mtext("(c)", line = -3, at = min(Ds), adj = 0, cex = 2.25)
legend(
    legendloc, bty = "n", lwd = 4, lty = 1, cex = legendlabelsize,
    col = "#33d17a", legend = "Sentinel nodes, n = 3"
)
## d: SA, n = 4
plot_ns(Ds, Y, font.lab = 3)
lines(Ds, y, lwd = 6, col = "#000000")
matlines(Ds, Y[, solns[[4]]$vs], lty = 1, lwd = 4, col = "#3584e4")
lines(Ds, rowMeans(Y[, solns[[4]]$vs]), lwd = 6, col = "#9141ac")
mtext("(d)", line = -3, at = min(Ds), adj = 0, cex = 2.25)
legend(
    legendloc, bty = "n", lwd = 4, lty = 1, cex = legendlabelsize,
    col = "#3584e4", legend = "Sentinel nodes, n = 4"
)
## e: GBB, including original and Chebyshev order 2 and 3
plot_ns(Ds, Y, font.lab = 3)
lines(Ds, y, lwd = 6, col = "#000000")
lines(Ds, GBB.obs, lwd = 6, col = "#986a44")
lines(Ds, GBB, lwd = 6, col = "#e01b24")
mtext("(e)", line = -3, at = min(Ds), adj = 0, cex = 2.25)
legend(
    legendloc, bty = "n", cex = legendlabelsize,
    lty = 1, lwd = 6,
    col = c("#986a44", "#e01b24"),
    legend = c(
        "GBB observable",
        "GBB approximation"
    )
)
## f: DART, including original and DART dim 2 and 3
plot_ns(Ds, Y, font.lab = 3)
lines(Ds, y, lwd = 6, col = "#000000")
lines(Ds, DART.obs, lwd = 6, col = "#e5a50a")
lines(Ds, DART, lwd = 6, col = "#ff007f")
mtext("(f)", line = -3, at = min(Ds), adj = 0, cex = 2.25)
legend(
    legendloc, bty = "n", cex = legendlabelsize,
    lty = 1, lwd = 6,
    col = c("#e5a50a", rep("#ff007f", 3)),
    legend = c(
        "DART observable",#, 1-D",
        "DART approximation"
    )
)
if(save_plots) dev.off()

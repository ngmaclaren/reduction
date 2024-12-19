library(optNS)

load("dolphin-demo-weights.RData")

unweighted <- new.env()
load("dolphin-demo.RData", envir = unweighted)

plot_ns <- function(Ds, Y, ...) {
    matplot(
        Ds, Y, type = "l", lty = 1, lwd = 0.5, col = "#babdb6",
        xlab = "D", ylab = "x", cex.lab = labelsize, cex.axis = labelsize, ...
    )
}

save_plots <- TRUE # FALSE
imgfile <- "../img/dolphin-demo-weights-v3.pdf"
labelsize <- 2.5
legendlabelsize <- 0.6*labelsize
legendloc <- "bottomright"

                                        # opt
get_error(solns)

Ds <- sdn::.doublewell$Ds

ht <- 13*2/3
wd <- 13
if(save_plots) {
    pdf(imgfile, height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
par(mfcol = c(2, 3), mar = c(5, 5, 1, 1), pty = "s")
## a: SA, n = 1
plot_ns(Ds, Y, font.lab = 3)
lines(Ds, y, lwd = 6, col = "#000000")
lines(Ds, with(unweighted, Y[, solns[[1]]$vs]), lwd = 6, col = "#9141ac")
lines(Ds, Y[, solns[[1]]$vs], lwd = 6, col = "#986a44") # 
lines(Ds, Y[, solns[[1]]$vs], lwd = 4, col = "#ff7800") 
mtext("(a)", line = -3, at = min(Ds), adj = 0, cex = 2.25)
legend(
    legendloc, bty = "n", lwd = c(2, 6, 6, 4), lty = 1, cex = legendlabelsize,
    col = c("#babdb6", "#000000", "#986a44", "#9141ac", "#ff7800"),
    ##legend = c("Node state", "Average state", "Sentinel node approximation", "Sentinel nodes, n = 1")
    legend = c("Node state", "Average state", "Weight-optimized", "Not weight-optimized", "Sentinel nodes, n = 1")
)
## c: SA, n = 3
plot_ns(Ds, Y, font.lab = 3)
lines(Ds, y, lwd = 6, col = "#000000")
lines(Ds, with(unweighted, rowMeans(Y[, solns[[3]]$vs])), lwd = 6, col = "#9141ac")
matlines(Ds, Y[, solns[[3]]$vs], lty = 1, lwd = 4, col = "#33d17a")
lines(Ds, apply(Y[, solns[[3]]$vs], 1, weighted.mean, w = solns[[3]]$ws), lwd = 6, col = "#986a44") # 9141ac
mtext("(c)", line = -3, at = min(Ds), adj = 0, cex = 2.25)
legend(
    legendloc, bty = "n", lwd = 4, lty = 1, cex = legendlabelsize,
    col = "#33d17a", legend = "Sentinel nodes, n = 3"
)
## b: SA, n = 2
plot_ns(Ds, Y, font.lab = 3)
lines(Ds, y, lwd = 6, col = "#000000")
lines(Ds, with(unweighted, rowMeans(Y[, solns[[2]]$vs])), lwd = 6, col = "#9141ac")
matlines(Ds, Y[, solns[[2]]$vs], lty = 1, lwd = 4, col = "#f6d32d")
lines(Ds, apply(Y[, solns[[2]]$vs], 1, weighted.mean, w = solns[[2]]$ws), lwd = 6, col = "#986a44") # 9141ac
mtext("(b)", line = -3, at = min(Ds), adj = 0, cex = 2.25)
legend(
    legendloc, bty = "n", lwd = 4, lty = 1, cex = legendlabelsize,
    col = "#f6d32d", legend = "Sentinel nodes, n = 2"
)
## d: SA, n = 4
plot_ns(Ds, Y, font.lab = 3)
lines(Ds, y, lwd = 6, col = "#000000")
lines(Ds, with(unweighted, rowMeans(Y[, solns[[4]]$vs])), lwd = 6, col = "#9141ac")
matlines(Ds, Y[, solns[[4]]$vs], lty = 1, lwd = 4, col = "#3584e4")
lines(Ds, apply(Y[, solns[[4]]$vs], 1, weighted.mean, w = solns[[4]]$ws), lwd = 6, col = "#986a44") # 9141ac
mtext("(d)", line = -3, at = min(Ds), adj = 0, cex = 2.25)
legend(
    legendloc, bty = "n", lwd = 4, lty = 1, cex = legendlabelsize,
    col = "#3584e4", legend = "Sentinel nodes, n = 4"
)
if(save_plots) dev.off()

load("dolphin-demo.RData")

save_plots <- TRUE # FALSE

library(latex2exp)
placelabel <- function(label, x, y, ...) {
    xlim <- par("usr")[1:2]
    ylim <- par("usr")[3:4]
    xpos <- xlim[1] + x*(xlim[2] - xlim[1])
    ypos <- ylim[1] + y*(ylim[2] - ylim[1])
    if(par("xlog")) xpos <- 10^xpos
    if(par("ylog")) ypos <- 10^ypos
    text(xpos, ypos, label, ...)
}

imgfile <- "../img/dolphin-demo.pdf"
colors <- list(
    nodestates = "#babdb6",
    systemstate = "#000000", 
    approximation = "#729fcf", 
    selectednodes = "#73d216", 
    GBB = "#e01b24", 
    DART = "#ff007f"
)
legendtext <- c("Node state", "Average state", "SN approximation", "Sentinel nodes", "GBB", "DART")

ht <- 12
wd <- 9
labelsize <- 2

if(save_plots) {
    pdf("../img/dolphin-demo.pdf", height = 12, width = wd)
} else {
    dev.new(height = ht, width = wd)
}

par(mfrow = c(3, 2), mar = c(5, 5, 1, 1))

plot_ns <- function(Ds, Y, color, labelsize) {
    matplot(
        Ds, Y, type = "l", lty = 1, lwd = 0.5, col = color,
        xlab = TeX(r"(D)", italic = TRUE), ylab = TeX(r"($x^*$)", italic = TRUE),
        cex.lab = labelsize, cex.axis = labelsize
    )
}
    
Ds <- sdn::.doublewell$Ds
for(i in seq_along(solns)) {
    plot_ns(Ds, Y, colors$nodestates, labelsize)

    lines(Ds, y, lty = 1, lwd = 6, col = colors$systemstate)

    if(ns[i] == 1) {
        lines(Ds, Y[, solns[[i]]$vs], lty = 1, lwd = 8, col = colors$approximation)
    } else {
        matlines(Ds, Y[, solns[[i]]$vs], lty = 1, lwd = 4, col = colors$selectednodes)
        lines(Ds, rowMeans(Y[, solns[[i]]$vs]), lty = 1, lwd = 8, col = colors$approximation)
    }

    placelabel(paste0("(", letters[i], ")"), 0.01, 0.99, adj = c(0, 1), cex = labelsize)

    if(i == 1) legend(
                   "bottomright", bty = "n", cex = 0.75*labelsize, lwd = 3, lty = 1,
                   col = unlist(colors), legend = legendtext
               )
}

## GBB
plot_ns(Ds, Y, colors$nodestates, labelsize)
lines(Ds, GBB.obs, lty = 1, lwd = 6, col = "#c17d11")
lines(Ds, y, lty = 1, lwd = 6, col = colors$systemstate)
lines(Ds, GBB[, 1], lty = 1, lwd = 8, col = colors$GBB)
placelabel(paste0("(", letters[5], ")"), 0.01, 0.99, adj = c(0, 1), cex = labelsize)
## DART
plot_ns(Ds, Y, colors$nodestates, labelsize)
lines(Ds, DART.obs, lty = 1, lwd = 6, col = "#c17d11")
lines(Ds, y, lty = 1, lwd = 6, col = colors$systemstate)
lines(Ds, DART[, 1], lty = 1, lwd = 8, col = colors$DART)
placelabel(paste0("(", letters[6], ")"), 0.01, 0.99, adj = c(0, 1), cex = labelsize)

if(save_plots) dev.off()

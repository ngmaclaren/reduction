load("dolphin-demo-weights.RData")

unweighted <- new.env()
load("dolphin-demo.RData", envir = unweighted)

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

imgfile <- "../img/dolphin-demo-weights.pdf"
colors <- list(
    nodestates = "#babdb6",
    systemstate = "#000000",
    weighted = "#e9b96e",
    approximation = "#3465a4",
    selectednodes = "#73d216"
)
legendtext <- c(
    "Node states",
    "System state",
    "Weight-optimized",
    "Not weight-optimized",
    "Selected nodes"
)

ht <- 8
wd <- 9
labelsize <- 2
palette("Tableau 10")

if(save_plots) {
    pdf(imgfile, height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}

par(mfrow = c(2, 2), mar = c(5, 5, 1, 1))

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
        lines(Ds, with(unweighted, Y[, solns[[1]]$vs]), lty = 1, lwd = 8,
              col = adjustcolor(colors$approximation, 1))

        lines(Ds, Y[, solns[[i]]$vs], lty = 1, lwd = 8, col = colors$weighted)
        
    } else {
        matlines(Ds, Y[, solns[[i]]$vs], lty = 1, lwd = 4, col = colors$selectednodes)
        lines(Ds, with(unweighted, rowMeans(Y[, solns[[i]]$vs])), lty = 1, lwd = 8,
              col = adjustcolor(colors$approximation, 1))

        lines(
            Ds,
            apply(Y[, solns[[i]]$vs], 1, weighted.mean, w = solns[[i]]$ws),
            lty = 1, lwd = 8, col = colors$weighted
        )

    }

    placelabel(paste0("(", letters[i], ")"), 0.01, 0.99, adj = c(0, 1), cex = labelsize)

    if(i == 1) legend(
                   "bottomright", bty = "n", cex = 0.75*labelsize, lwd = 3, lty = 1,
                   col = unlist(colors), legend = legendtext
               )
}

if(save_plots) dev.off()

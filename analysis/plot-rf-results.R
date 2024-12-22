library(sfsmisc)
palette("Okabe-Ito")
red <- palette.colors(2, "R4")[2]

df <- read.csv("../data/rf-importances.csv")
df <- df[order(df$N), ]
xcol <- c("k", "knn", "lcl", "cc", "bc", "kcore")

with(list(df = subset(df, network != "er")), {
    dev.new(height = 7, width = 8)
    par(mar = c(5, 5, 1, 5), pty = "s")
    matplot(
        df$N, df[, xcol], type = "o", log = "x",
        pch = 21, lty = 1, lwd = 2, col = seq(length(xcol))+1, bg = "white", cex = 1.5,
        xlab = "", ylab = "", cex.lab = 2, axes = FALSE
    )
    par(new = TRUE)
    plot(
        df$N, df$r2, type = "l", col = red, lwd = 6, ylim = c(-1, 1), log = "x",
        xlab = "", ylab = "", axes = FALSE
    )
    eaxis(1, at = axTicks(1)[c(FALSE, TRUE)], cex.axis = 1.75)
    eaxis(2, cex.axis = 1.75)
    eaxis(
        4, cex.axis = 1.75, col.axis = red,
        col = red, small.args = list(col = red)
    )
    mtext("Coefficient of Determination", side = 4, line = 4, col = red, cex = 1.75)
    box()
    title(xlab = "N", font.lab = 3, cex.lab = 2)
    title(ylab = "Importance", cex.lab = 2, line = 4)
    legend(
        "topright", bty = "n",
        pch = 21, pt.cex = 1.5, lty = 1, lwd = 2, col = seq(length(xcol))+1, pt.bg = "white",
        legend = c(
            "Degree", "Avg. nearest neighbor deg.", "Local clustering", "Closeness", "Betweenness", "K-core"
        )
    )
})

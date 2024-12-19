library(sfsmisc)
setwd("~/Documents/reduction/shell/output/")
save_plots <- TRUE # FALSE

filenames <- list.files(pattern = ".csv")
csvs <- lapply(filenames, read.csv)

df <- do.call(rbind, csvs)
df$dynamics <- factor(df$dynamics, levels = c("doublewell", "mutualistic", "SIS", "genereg"))

agg <- aggregate(runtime ~ network + N + dynamics, data = df, FUN = max) # mean
agg$color <- as.integer(agg$dynamics)
agg$shape <- as.integer(agg$dynamics)

guide <- lm(log(runtime) ~ log(N), data = agg)

if(save_plots) {
    pdf("~/Documents/reduction/img/SA-runtime.pdf")
} else {
    dev.new()
}
par(mar = c(5, 5, 1, 1))
palette("Set 1")
plot(
    runtime ~ N, data = agg, type = "p",
    col = agg$color, cex = 2, lwd = 2, pch = agg$shape,
    xlab = "Number of nodes", ylab = "Wall time",
    cex.lab = 1.75, log = "xy",
    axes = FALSE
)
with(list(x = seq(500, 20000)), lines(x, exp(coef(guide)[1])*x^2)) # force x^2
text(200, 100, labels = expression(. %prop% N^2), cex = 1.5, adj = 0)
box()
eaxis(1, at = axTicks(1)[c(FALSE, TRUE)], cex.axis = 1.5)
eaxis(2, cex.axis = 1.5)
legend(
    "topleft", bty = "n",
    legend = c("Double-well", "Mutualistic species", "SIS", "Gene-regulatory"),
    pch = 1:4, col = 1:4, pt.cex = 1.5, pt.lwd = 2, cex = 1.25
)
if(save_plots) dev.off()

## m <- lm(runtime ~ poly(N, 2, raw = TRUE), data = agg)

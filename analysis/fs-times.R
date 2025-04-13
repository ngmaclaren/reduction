library(sfsmisc)

save_plots <- TRUE # FALSE

df <- read.csv("../data/fs-times-20250125.csv")

rdf <- reshape(
    data = df,
    varying = c("dw", "sis", "gr", "ms"), v.names = "runtime", timevar = "dynamics",
    times = c("doublewell", "SIS", "genereg", "mutualistic"),
    new.row.names = seq(1000), direction = "long"
)

rdf$dynamics <- factor(rdf$dynamics, levels = c("doublewell", "mutualistic", "SIS", "genereg"))

rdf$color <- as.integer(rdf$dynamics)
rdf$shape <- as.integer(rdf$dynamics)

guides <- lapply(levels(rdf$dynamics), function(dyn) {
    lm(log(runtime) ~ log(N), data = rdf, subset = dynamics == dyn)
})

if(save_plots) {
    pdf("~/Documents/reduction/img/FS-runtime.pdf")
} else {
    dev.new()
}
par(mar = c(5, 5, 1, 1))
palette("Set 1")
plot(
    runtime ~ N, data = rdf, type = "p",
    col = rdf$color, cex = 2, lwd = 2, pch = rdf$shape,
    xlab = "Number of nodes", ylab = "Wall time",
    cex.lab = 1.75, log = "xy",
    axes = FALSE
)
with(list(x = seq(500, 20000)), lines(x, exp(coef(guides[[1]])[1])*x^2, col = "black", lwd = 1)) # force x^2
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

## m <- lm(runtime ~ poly(N, 2, raw = TRUE), data = rdf)

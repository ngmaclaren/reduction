library(parallel)
ncores <- detectCores() - 1
library(igraph)
library(deSolve)
source("../src/functions.R")

network <- "dolphin"

load("../data/dolphin.rda")
g <- get(network)
N <- vcount(g)
A <- as_adj(g, "both", sparse = FALSE)
k <- degree(g)

x.init <- 1
B <- 1
f <- 1
h <- 2
                                        # This is the change
Ds <- seq(1, 0.001, length.out = 1000)

Y <- solve_genereg(x = rep(x.init, N), B = B, f = f, h = h, Ds = Ds, A = A)

                                        # color the nodes? 
colors <- hcl.colors(10)
colorramp <- colorRamp(colors)
nc <- colorramp(k/max(k))
nodecolors <- apply(nc, 1, function(r) rgb(r[1], r[2], r[3], maxColorValue = 255))
colorbar <- function(colors, rng, title = "") {
    ##crp <- colorRampPalette(colors)
    z <- matrix(1:100, nrow = 1)
    x <- 1
    y <- seq(rng[1], rng[2], length.out = 100)
    par(mar = c(6, 0, 3, 3.5) + .5)
    image(x, y, z, col = colors,#crp(100),
          axes = FALSE, xlab = "", ylab = "")
    axis(4, at = pretty(rng), labels = pretty(rng), cex.axis = ticksize)
    mtext(title, 4, line = 2.5, cex = labelsize)
    box()
}


pdf("genereg-bifurcation-diagram.pdf", width = 8)
labelsize <- 1.75
ticksize <- 1.75
nf <- layout(matrix(c(1, 2), nrow = 1, ncol = 2), c(4.1, 0.9), 5)
matplot(
    Ds, Y, type = "l", lwd = .75, lty = 1, col = nodecolors,
    xlab = "D", ylab = "x", cex.axis = ticksize, cex.lab = labelsize
)
colorbar(hcl.colors(100), range(k), "Degree")
dev.off()

                                        # Focus
Ds <- seq(0.35, 0.25, length.out = 1000)
Y <- solve_genereg(x = rep(x.init, N), B = B, f = f, h = h, Ds = Ds, A = A)
pdf("genereg-bifurcation-diagram-focused.pdf", width = 8)
labelsize <- 1.75
ticksize <- 1.75
nf <- layout(matrix(c(1, 2), nrow = 1, ncol = 2), c(4.1, 0.9), 5)
matplot(
    Ds, Y, type = "l", lwd = .75, lty = 1, col = nodecolors,
    xlab = "D", ylab = "x", cex.axis = ticksize, cex.lab = labelsize
)
colorbar(hcl.colors(100), range(k), "Degree")
dev.off()

                                        # Color by community
part <- cluster_edge_betweenness(g)
pdf("dolphin-communities.pdf")
plot(part, g, palette = "Tableau 10")
dev.off()
pdf("genereg-bifurcation-diagram-focused-communities.pdf")
palette("Tableau 10")
matplot(
    Ds, Y, type = "l", lwd = 1, lty = 1, col = membership(part),
    xlab = "D", ylab = "x", cex.axis = ticksize, cex.lab = labelsize
)
legend("topleft", legend = sort(unique(membership(part))), col = sort(unique(membership(part))),
       lty = 1, lwd = 3, bty = "n",
       title = "Community")
dev.off()

                                        # ?!
                                        # Do the same for SIS
sis <- new.env()
with(sis, {
    x.init <- SIS_parms$x.init
    mu <- 1
    Ds <- SIS_parms$Ds
    Y <- solve_SIS(rep(x.init, N), mu, Ds, A)
    pdf("SIS-bifurcation-diagram.pdf", width = 8)
    nf <- layout(matrix(c(1, 2), nrow = 1, ncol = 2), c(4.1, 0.9), 5)
    matplot(
        Ds, Y, type = "l", lwd = .75, lty = 1, col = nodecolors,
        xlab = "D", ylab = "x", cex.axis = ticksize, cex.lab = labelsize
    )
    colorbar(hcl.colors(100), range(k), "Degree")
    dev.off()
})

with(sis, {
    Ds <- seq(0.1, 0.3, length.out = 1000)
    Y <- solve_SIS(rep(x.init, N), mu, Ds, A)
    pdf("SIS-bifurcation-diagram-focused.pdf", width = 8)
    nf <- layout(matrix(c(1, 2), nrow = 1, ncol = 2), c(4.1, 0.9), 5)
    matplot(
        Ds, Y, type = "l", lwd = .75, lty = 1, col = nodecolors,
        xlab = "D", ylab = "x", cex.axis = ticksize, cex.lab = labelsize
    )
    colorbar(hcl.colors(100), range(k), "Degree")
    dev.off()

    pdf("SIS-bifurcation-diagram-focused-communities.pdf")
    palette("Tableau 10")
    matplot(
        Ds, Y, type = "l", lwd = 1, lty = 1, col = membership(part),
        xlab = "D", ylab = "x", cex.axis = ticksize, cex.lab = labelsize
    )
    legend("topleft", legend = sort(unique(membership(part))), col = sort(unique(membership(part))),
           lty = 1, lwd = 3, bty = "n",
           title = "Community")
    dev.off()
})

                                        # ??!!
                                        # Now try doublewell
dw <- new.env()
with(dw, {
    x.init <- 1
    r <- c(1, 3, 5)
    Ds <- seq(0.001, 1, length.out = 1000)
    Y <- solve_doublewell(rep(x.init, N), r, Ds, A)
    pdf("doublewell-bifurcation-diagram.pdf", width = 8)
    nf <- layout(matrix(c(1, 2), nrow = 1, ncol = 2), c(4.1, 0.9), 5)
    matplot(
        Ds, Y, type = "l", lwd = .75, lty = 1, col = nodecolors,
        xlab = "D", ylab = "x", cex.axis = ticksize, cex.lab = labelsize
    )
    colorbar(hcl.colors(100), range(k), "Degree")
    dev.off()

    pdf("doublewell-bifurcation-diagram-communities.pdf")
    palette("Tableau 10")
    matplot(
        Ds, Y, type = "l", lwd = 1, lty = 1, col = membership(part),
        xlab = "D", ylab = "x", cex.axis = ticksize, cex.lab = labelsize
    )
    legend("topleft", legend = sort(unique(membership(part))), col = sort(unique(membership(part))),
           lty = 1, lwd = 3, bty = "n",
           title = "Community")
    dev.off()
})

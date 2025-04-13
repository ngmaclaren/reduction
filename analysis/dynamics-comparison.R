library(sfsmisc)
library(optNS)
library(igraph)

save_plots <- TRUE # FALSE

imgfile <- "../img/dynamics-comparison-v2.pdf"

Ylist <- readRDS("../data/fullstate-dolphin.rds") # still fixing dolphin network

ns.types <- c("opt", "fixed", "rand")
legendtext <- c("Optimized", "Degree-preserving", "Random")

make_plotdata <- function(dynamicsA, dynamicsB, ns.types, Ylist) {
    A <- readRDS(paste0("../data/ns-dolphin_", dynamicsA, ".rds"))[ns.types]
    A.Y <- Ylist[[dynamicsA]]
    A.y <- rowMeans(A.Y)

    B <- readRDS(paste0("../data/ns-dolphin_", dynamicsB, ".rds"))[ns.types]
    B.Y <- Ylist[[dynamicsB]]
    B.y <- rowMeans(B.Y)

    list(
        x = lapply(A, get_error),
        y = lapply(A, function(ns) sapply(ns, function(x) obj_fn(x$vs, B.y, B.Y)))
    )
}

dyns <- c("doublewell", "mutualistic", "SIS", "genereg")
conds <- expand.grid(dyns, dyns, stringsAsFactors = FALSE)
colnames(conds) <- c("dynamicsB", "dynamicsA")
labnames <- c("Double-well", "Mutualistic species", "SIS", "Gene-regulatory")
labnames <- expand.grid(labnames, labnames, stringsAsFactors = FALSE)
colnames(labnames) <- c("nameB", "nameA")
conds <- cbind(conds, labnames)
conds <- conds[conds$dynamicsA != conds$dynamicsB, ]
conds <- conds[, c("dynamicsA", "dynamicsB", "nameA", "nameB")]


show_axis_labels <- FALSE # TRUE
unit <- 2.5
ht <- 4*unit
wd <- 4*unit
labelsize <- 1.25
colors <- list(opt = "#3584e4", fixed = "#ff7800", rand = "#33d17a")
pchs <- list(opt = 1, fixed = 2, rand = 0)

if(save_plots) {
    pdf(imgfile, height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}

par(mfcol = c(4, 4))
if(show_axis_labels) {
    par(mar = c(4, 4.1, 1, 0.9))
} else {
    par(mar = c(2.75, 2.75, 0.25, 0.25))
}

eqs <- list(
    eqno = 1:4,
    dynamics = dyns,
    name = c("Double-well", "Mutualistic species", "SIS", "Gene-regulatory"),
    model = c(
        expression(frac(d*x[i], d*t) == -(x[i] - r[1])*(x[i] - r[2])*(x[i] - r[3]) + D*sum(a[ij]*x, j=1,N)),
        expression(frac(d*x[i], d*t) == B + x[i]*bgroup("(", 1 - frac(x[i], K), ")")*bgroup("(", frac(x[i], C) - 1, ")") + D*sum(a[ij]*frac(x[i]*x[j], tilde(D) + E*x[i] + H*x[j]), j=1, N)),
        expression(frac(d*x[i], d*t) == -mu*x[i] + lambda*sum(a[ij]*(1 - x[i])*x[j], j=1, N)),
        expression(frac(d*x[i], d*t) == -B*x[i]^f + D*sum(a[ij]*frac(x[j]^h, 1 + x[j]^h), j=1, N))
    )
)

i <- 1
eq <- 1
for(pltno in seq(16)) {
    if(any(pltno == 1 + seq(0, 15, by = 5))) {
        plot(NULL, axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1))
        ## next
        text(0.5, 0.66, labels = eqs$name[eq], cex = 1.5, font = 2)
        text(0.5, 0.33, labels = eqs$model[eq], cex = switch(eq, 1, 0.75, 1, 1))
        eq <- eq + 1
        next
    }
    
    plotdata <- make_plotdata(conds[i, "dynamicsA"], conds[i, "dynamicsB"], ns.types, Ylist)
    x <- plotdata$x
    y <- plotdata$y
    
    xlim <- range(unlist(x))
    ylim <- range(unlist(y))

    plot(NULL, xlim = xlim, ylim = ylim, xlab = "", ylab = "", xaxt = "n", yaxt = "n", log = "xy")
    eaxis(1, n.axp = 1, cex.axis = labelsize)
    eaxis(2, n.axp = 1, cex.axis = labelsize)

    if(show_axis_labels) {
        title(xlab = conds[i, "nameA"], cex.lab = labelsize)
        title(ylab = conds[i, "nameB"], cex.lab = labelsize)
    }
    for(nst in ns.types) points(x[[nst]], y[[nst]], col = adjustcolor(colors[[nst]], .25), pch = pchs[[nst]])
    for(nst in ns.types) points(mean(x[[nst]]), mean(y[[nst]]), col = colors[[nst]],
                                pch = pchs[[nst]], lwd = 3, cex = 3)
   
    if(i == 1) {
        legend(
            "bottomright", bty = "n", col = unlist(colors), pch = unlist(pchs), pt.lwd = 2, pt.cex = 2,
            cex = labelsize, legend = legendtext
        )
    }

    i <- i + 1
}

if(save_plots) dev.off()

## library(latex2exp)
ns <- readRDS("../data/ns-dolphin_doublewell.rds")$opt # SIS
S <- ns[[which.min(get_error(ns))]]
g <- readRDS("../data/dolphin.rds")
## print(V(g)$name[S$vs])
print(S$vs)
## V(g)$color <- NA # V(g)[degree(g) == 11]$color <- 3
## V(g)[S$vs]$color <- 1
## dev.new(height = 10, width = 10)
## plot(g, vertex.size = 8)


Ds <- sdn::.doublewell$Ds
ht2 <- 5
wd2 <- 5
ls2 <- 2

xicolor <- "#c0bfbc"
vscolor <- "#3584e4"
appcolor <- "#9141ac"

                                        # on SIS
if(save_plots) {
    pdf("../img/dynamics-comparison-SISbifplot-v2.pdf", height = ht2, width = wd2)
} else {
    dev.new(height = ht2, width = wd2)
}
par(mar = c(5, 5, 1, 1))
matplot(
    Ds, Ylist$SIS, type = "l", lty = 1, lwd = 0.5, col = xicolor, ylim = c(0, 1),
    ## xlab = TeX(r"(D)", italic = TRUE), ylab = TeX(r"($x^*$)", italic = TRUE),
    xlab = "D", ylab = "x*",
    font.lab = 3, cex.lab = ls2, cex.axis = ls2
)
lines(Ds, rowMeans(Ylist$SIS), lty = 1, lwd = 6, col = "#000000")
matlines(Ds, Ylist$SIS[, S$vs], lty = 1, lwd = 4, col = vscolor)
lines(Ds, rowMeans(Ylist$SIS[, S$vs]), lty = 1, lwd = 8, col = appcolor)
if(save_plots) dev.off()

                                        # on Double-well
if(save_plots) {
    pdf("../img/dynamics-comparison-doublewellbifplot-v2.pdf", height = ht2, width = wd2)
} else {
    dev.new(height = ht2, width = wd2)
}
par(mar = c(5, 5, 1, 1))
matplot(
    Ds, Ylist$doublewell, type = "l", lty = 1, lwd = 0.5, col = xicolor,
    ## xlab = TeX(r"(D)", italic = TRUE), ylab = TeX(r"($x^*$)", italic = TRUE),
    xlab = "D", ylab = "x*",
    font.lab = 3, cex.lab = ls2, cex.axis = ls2
)
lines(Ds, rowMeans(Ylist$doublewell), lty = 1, lwd = 6, col = "#000000")
matlines(Ds, Ylist$doublewell[, S$vs], lty = 1, lwd = 4, col = vscolor)
lines(Ds, rowMeans(Ylist$doublewell[, S$vs]), lty = 1, lwd = 8, col = appcolor)
if(save_plots) dev.off()

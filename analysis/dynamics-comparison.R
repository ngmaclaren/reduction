library(sfsmisc)
library(optNS)

save_plots <- TRUE # FALSE

imgfile <- "../img/dynamics-comparison.pdf"

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

unit <- 2.5
ht <- 4*unit
wd <- 3*unit
labelsize <- 1.5

colors <- list(opt = "#3584e4", fixed = "#ff7800", rand = "#33d17a")
pchs <- list(opt = 1, fixed = 2, rand = 0)

dyns <- c("doublewell", "mutualistic", "genereg", "SIS")
conds <- expand.grid(dyns, dyns, stringsAsFactors = FALSE)
colnames(conds) <- c("dynamicsB", "dynamicsA")
labnames <- c("Double-well", "Mutualistic species", "Gene-regulatory", "SIS")
labnames <- expand.grid(labnames, labnames, stringsAsFactors = FALSE)
colnames(labnames) <- c("nameB", "nameA")
conds <- cbind(conds, labnames)
conds <- conds[conds$dynamicsA != conds$dynamicsB, ]
conds <- conds[, c("dynamicsA", "dynamicsB", "nameA", "nameB")]


if(save_plots) {
    pdf(imgfile, height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}

par(mfrow = c(4, 3), mar = c(4, 4.1, 1, 0.9))

for(i in seq(nrow(conds))) {
    plotdata <- make_plotdata(conds[i, "dynamicsA"], conds[i, "dynamicsB"], ns.types, Ylist)
    x <- plotdata$x
    y <- plotdata$y
    
    xlim <- range(unlist(x))
    ylim <- range(unlist(y))

    plot(NULL, xlim = xlim, ylim = ylim, xlab = "", ylab = "", xaxt = "n", yaxt = "n", log = "xy")
    eaxis(1, n.axp = 1, cex.axis = labelsize)
    eaxis(2, n.axp = 1, cex.axis = labelsize)
    title(xlab = conds[i, "nameA"], cex.lab = labelsize)
    title(ylab = conds[i, "nameB"], cex.lab = labelsize)
    for(nst in ns.types) points(x[[nst]], y[[nst]], col = adjustcolor(colors[[nst]], .25), pch = pchs[[nst]])
    for(nst in ns.types) points(mean(x[[nst]]), mean(y[[nst]]), col = colors[[nst]],
                                pch = pchs[[nst]], lwd = 3, cex = 3)
   
    if(i == 1) {
        legend(
            "bottomright", bty = "n", col = unlist(colors), pch = unlist(pchs), pt.lwd = 2, pt.cex = 2,
            cex = labelsize, legend = legendtext
        )
    }
}

if(save_plots) dev.off()

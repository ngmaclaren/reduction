library(sfsmisc)
library(optNS)

save_plots <- TRUE # FALSE
useall <- "no" # "yes" # broken right now. Need to add knn-constrained
useweights <- "no" # "yes"

A.dynamics <- "doublewell"
B.dynamics <- "SIS"

imgfile <- switch(useweights, no = "../img/dynamics-comparison.pdf", yes = "../img/dynamics-comparison_w.pdf")

Ylist <- readRDS("../data/fullstate-dolphin.rds")

A <- readRDS(switch(useweights, no = "../data/ns-dolphin_doublewell.rds", yes = "../data/ns-dolphin_doublewell_w.rds"))
A.Y <- Ylist$doublewell
A.y <- rowMeans(A.Y)
B <- readRDS(switch(useweights, no = "../data/ns-dolphin_SIS.rds", yes = "../data/ns-dolphin_SIS_w.rds"))
B.Y <- Ylist$SIS
B.y <- rowMeans(B.Y)

weightsflag <- switch(useweights, no = FALSE, yes = TRUE)
x1 <- lapply(A, get_error)
y1 <- lapply(A, function(ns) sapply(ns, function(x) obj_fn(x$vs, B.y, B.Y, weightsflag, ws = x$ws)))

x2 <- lapply(B, get_error)
y2 <- lapply(B, function(ns) sapply(ns, function(x) obj_fn(x$vs, A.y, A.Y, weightsflag, ws = x$ws)))

usethese <- switch(useall, no = 1:3, yes = seq_along(A))
legendtext <- c("Optimized", "Degree-preserving", "Random", "Constrained", "Quantiled", "Community-based")

ht <- 7
wd <- 14
palette("Set 1")
## palette("Tableau 10")
labelsize <- 1.5

if(save_plots) {
    pdf(imgfile, height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}

par(mfrow = c(1, 2), mar = c(4, 4.1, 1, 0.9))

colors <- c("#3584e4", "#ff7800", "#33d17a")
pchs <- c(1, 2, 0)

xlim1 <- range(unlist(x1[usethese]))
ylim1 <- range(unlist(y1[usethese]))
plot(NULL, xlim = xlim1, ylim = ylim1, xlab = "", ylab = "", xaxt = "n", yaxt = "n", log = "xy")
eaxis(1, n.axp = 1, cex.axis = labelsize)
eaxis(2, n.axp = 1, cex.axis = labelsize)
title(xlab = "Double-well", cex.lab = labelsize)
title(ylab = "SIS", cex.lab = labelsize)
for(i in usethese) {
    points(x1[[i]], y1[[i]], col = colors[i], pch = pchs[i])
    points(mean(x1[[i]]), mean(y1[[i]]), col = colors[i], pch = pchs[i], lwd = 3, cex = 3)
}
legend(
    "bottomright", bty = "n", col = colors, pch = pchs, pt.lwd = 2, pt.cex = 2, cex = labelsize,
    legend = legendtext[usethese]
)

xlim2 <- range(unlist(x2[usethese]))
ylim2 <- range(unlist(y2[usethese]))
plot(NULL, xlim = xlim2, ylim = ylim2, xlab = "", ylab = "", xaxt = "n", yaxt = "n", log = "xy")
eaxis(1, n.axp = 1, cex.axis = labelsize)
eaxis(2, n.axp = 1, cex.axis = labelsize)
title(ylab = "Double-well", cex.lab = labelsize)
title(xlab = "SIS", cex.lab = labelsize)
for(i in usethese) {
    points(x2[[i]], y2[[i]], col = colors[i], pch = pchs[i])
    points(mean(x2[[i]]), mean(y2[[i]]), col = colors[i], pch = pchs[i], lwd = 3, cex = 3)
}

if(save_plots) dev.off()

## Pass in the optimized node weights from A, don't try to reoptimize and don't ignore

## The code below should be correct, but need to check somehow

library(sfsmisc)
library(optNS)

useall <- "no" # "yes"
useweights <- "yes" # "no"

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

dev.new(width = 14)
par(mfrow = c(1, 2), mar = c(4, 4.1, 1, 0.9))
palette("Set 1")
labelsize <- 1.5

xlim1 <- range(unlist(x1[usethese]))
ylim1 <- range(unlist(y1[usethese]))
plot(NULL, xlim = xlim1, ylim = ylim1, xlab = "Double-well", ylab = "SIS", log = "xy",
     cex.axis = labelsize, cex.lab = labelsize)
for(i in usethese) {
    points(x1[[i]], y1[[i]], col = i, pch = i)
    points(mean(x1[[i]]), mean(y1[[i]]), col = i, pch = i, lwd = 3, cex = 3)
}
legend(
    "bottomright", bty = "n", col = usethese, pch = usethese, pt.lwd = 2, pt.cex = 2, cex = labelsize,
    legend = legendtext[usethese]
)

xlim2 <- range(unlist(x2[usethese]))
ylim2 <- range(unlist(y2[usethese]))
plot(NULL, xlim = xlim2, ylim = ylim2, ylab = "Double-well", xlab = "SIS", log = "xy",
     cex.axis = labelsize, cex.lab = labelsize)
for(i in usethese) {
    points(x2[[i]], y2[[i]], col = i, pch = i)
    points(mean(x2[[i]]), mean(y2[[i]]), col = i, pch = i, lwd = 3, cex = 3)
}


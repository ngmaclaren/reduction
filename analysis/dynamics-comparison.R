library(sfsmisc)
library(optNS)

network <- "ba"
A.dynamics <- "doublewell"
B.dynamics <- "SIS"
useall <- "no" # "yes"
imgfile <- "../img/dynamics-comparison.pdf"

Ylist <- readRDS("../data/fullstate-ba.rds")

A <- readRDS("../data/ns-ba_doublewell.rds")
A.Y <- Ylist$doublewell
A.y <- rowMeans(A.Y)
B <- readRDS("../data/ns-ba_SIS.rds")
B.Y <- Ylist$SIS
B.y <- rowMeans(B.Y)

x1 <- lapply(A, get_error)
y1 <- lapply(A, function(ns) sapply(ns, function(x) obj_fn(x$vs, B.y, B.Y)))

x2 <- lapply(B, get_error)
y2 <- lapply(B, function(ns) sapply(ns, function(x) obj_fn(x$vs, A.y, A.Y)))

if(useall == "no") usethese <- 1:3 else usethese <- seq_along(A)

dev.new(width = 14)
par(mfrow = c(1, 2))
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

xlim2 <- range(unlist(x2[usethese]))
ylim2 <- range(unlist(y2[usethese]))
plot(NULL, xlim = xlim2, ylim = ylim2, ylab = "Double-well", xlab = "SIS", log = "xy",
     cex.axis = labelsize, cex.lab = labelsize)
for(i in usethese) {
    points(x2[[i]], y2[[i]], col = i, pch = i)
    points(mean(x2[[i]]), mean(y2[[i]]), col = i, pch = i, lwd = 3, cex = 3)
}


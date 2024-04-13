library(igraph)
library(optNS)
library(sfsmisc)
palette("Tableau 10")

save_plots <- FALSE # TRUE

dolphin <- readRDS("../data/dolphin.rds")
dolphin.N <- vcount(dolphin)
dolphin.k <- degree(dolphin)
dolphin.breaks <- seq(max(dolphin.k))
dolphin.hist <- hist(dolphin.k, breaks = dolphin.breaks, plot = FALSE)
dolphin.ns <- 1:4
dolphin.opts <- readRDS("../data/dolphin-1-to-12.rds")


                                        # lets try BA first.
ba <- readRDS("../data/ba.rds")
ba.N <- vcount(ba)
ba.k <- degree(ba)
ba.hist <- hist(log10(ba.k), plot = FALSE)
ba.hist$counts <- log10(ba.hist$counts)
ba.breaks <- ba.hist$breaks
ba.ns <- c(1, 2, 3, floor(log(ba.N)))
ba.opts <- readRDS("../data/ba-1-to-12.rds")

##plot(ba.hist)

prep_for_hist <- function(ns.list, logtransf = FALSE) {
    ks <- as.matrix(get_ks(ns.list))
    if(nrow(ks) < ncol(ks)) ks <- t(ks)
    ks <- as.data.frame(ks)
    colnames(ks) <- paste0("n", seq(ncol(ks)))
    dat <- reshape(ks, varying = colnames(ks), v.names = "ki", timevar = "n", times = seq(ncol(ks)),
                   new.row.names = 1:10000, direction = "long")
    if(logtransf) dat$logki <- log10(dat$ki)

    dat
}

dolphin.data <- lapply(dolphin.opts, prep_for_hist)
ba.data <- lapply(ba.opts, prep_for_hist, logtransf = TRUE)

make_histpanel <- function(plotdata, breaks, logtransf = FALSE) {
    ns <- unique(plotdata$n)
    if(logtransf) var <- "logki" else var <- "ki"
    hist(plotdata[, var], breaks = breaks, col = max(ns), main = "", xlab = "",
         ylab = switch(logtransf + 1, "", NULL),
         cex.lab = labelsize, axes = FALSE)
    box()
    ## xaxt = "n", yaxt = "n")
    if(logtransf) {
                                        # global!
        ##axis(1, at = axTicks(1)[c(FALSE, TRUE)], labels = 10^axTicks(1)[c(FALSE, TRUE)], cex.axis = labelsize)
        eaxis(1, cex.axis = labelsize)
        ## axis(2, at = axTicks(2), labels = 10^axTicks(1), cex.axis = labelsize)
        ##axis(2, cex.axis = labelsize)
        eaxis(2, cex.axis = labelsize)
    } else {
        eaxis(1, cex.axis = labelsize)
        eaxis(2, cex.axis = labelsize)
    }
    for(n in rev(ns[-length(ns)])) {
        hist(plotdata[, var][plotdata$n %in% seq(n)], breaks = breaks, col = n, add = TRUE)
    }
}

ht <- 10
wd <- 12
labelsize <- 2
if(save_plots) {
    pdf("../img/degree-sequences.pdf", height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
par(mfcol = c(5, 2), mar = c(4, 4, 2, 1) + 0.5)
## plot(dolphin.hist, main = "", xlab = "", cex.lab = labelsize, cex.axis = labelsize)
plot(dolphin.hist, main = "", xlab = "", axes = FALSE, cex.lab = labelsize)
box()
eaxis(1, cex.axis = labelsize)
eaxis(2, cex.axis = labelsize)
for(i in dolphin.ns) make_histpanel(dolphin.data[[i]], dolphin.breaks)
title(xlab = "Degree", cex.lab = labelsize)
## plot(ba.hist, main = "", xlab = "", cex.lab = labelsize, xaxt = "n", yaxt = "n")
## axis(1, at = axTicks(1)[c(FALSE, TRUE)], labels = 10^axTicks(1)[c(FALSE, TRUE)], cex.axis = labelsize)
## axis(2, at = axTicks(2), labels = 10^axTicks(2), cex.axis = labelsize)
plot(ba.hist, main = "", xlab = "", ylab = "", cex.lab = labelsize, axes = FALSE, log = "y")
box()
eaxis(1, cex.axis = labelsize)
eaxis(2, cex.axis = labelsize)
for(i in ba.ns) make_histpanel(ba.data[[i]], ba.breaks, logtransf = TRUE)
title(xlab = "Degree", cex.lab = labelsize)
if(save_plots) dev.off()

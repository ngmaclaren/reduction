library(igraph)
library(optNS)

save_plots <- TRUE # FALSE

dolphin <- readRDS("../data/dolphin.rds")
dolphin.N <- vcount(dolphin)
dolphin.k <- degree(dolphin)
dolphin.kposs <- seq(1, max(dolphin.k))

ba <- readRDS("../data/ba.rds")
ba.N <- vcount(ba)
ba.k <- degree(ba)
ba.kposs <- seq(1, max(ba.k)) # as.numeric(names(table(ba.k))) may work here as well

dolphin.opts <- readRDS("../data/dolphin-1-to-12.rds")
ba.opts <- readRDS("../data/ba-4.rds")

## need to make the barplot-ready data
prep_for_barplot <- function(nslist, kposs) { # a list of optimized node sets, possible degrees
    x <- as.matrix(get_ks(nslist))
    if(dim(x)[1] > dim(x)[2]) x <- t(x)
    y <- apply(x, 1, function(ki) sapply(kposs, function(k) sum(k == ki)))
    if(dim(y)[1] > dim(y)[2]) y <- t(y)
    return(y)
}

dolphin.data <- lapply(dolphin.opts[1:4], prep_for_barplot, dolphin.kposs)
ba.data <- lapply(ba.opts[1:4], prep_for_barplot, ba.kposs)

palette("Tableau 10")
ht <- 10
wd <- 12
labelsize <- 2
if(save_plots) {
    pdf("../img/degree-sequences.pdf", height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
par(mfcol = c(5, 2), mar = c(4, 4, 2, 0) + 0.5)
                                        # Dolphin
barplot(table(dolphin.k), names.arg = dolphin.kposs, col = "gray60",
        cex.axis = labelsize, cex.names = labelsize, cex.lab = labelsize,
        xlab = "", ylab = "Frequency")
for(i in 1:3) {
    barplot(dolphin.data[[i]], names.arg = dolphin.kposs, col = seq_len(nrow(dolphin.data[[i]])),
            cex.axis = labelsize, cex.names = labelsize, cex.lab = labelsize,
            xlab = "", ylab = "Frequency")
}
barplot(dolphin.data[[4]], names.arg = dolphin.kposs, col = seq_len(nrow(dolphin.data[[4]])),
        cex.axis = labelsize, cex.names = labelsize, cex.lab = labelsize,
        xlab = "Degree", ylab = "Frequency")
                                        # BA
barplot(sapply(ba.kposs, function(k) sum(ba.k == k)), names.arg = ba.kposs, col = "gray60",
        cex.axis = labelsize, cex.names = labelsize, cex.lab = labelsize,
        xlab = "", ylab = "")
for(i in 1:3) {
    barplot(ba.data[[i]], names.arg = ba.kposs, col = seq_len(nrow(ba.data[[i]])),
            cex.axis = labelsize, cex.names = labelsize, cex.lab = labelsize,
            xlab = "", ylab = "")
}
barplot(ba.data[[4]], names.arg = ba.kposs, col = seq_len(nrow(ba.data[[4]])),
        cex.axis = labelsize, cex.names = labelsize, cex.lab = labelsize,
        xlab = "Degree", ylab = "")

if(save_plots) dev.off()

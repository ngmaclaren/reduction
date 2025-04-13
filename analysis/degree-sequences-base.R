library(igraph)
library(optNS)

save_plots <- FALSE # TRUE

get_plims <- function(dims) {
    ## want to divide the unit square into pieces based on dim
    stopifnot(length(dims) == 2)
    stopifnot(all(floor(dims) == dims))
    list(horiz = seq(0, 1, by = 1/dims[2]), vert = rev(seq(0, 1, by = 1/dims[1])))
}
get_ploc <- function(loc, dims, byrow = FALSE) {
    lyt <- matrix(seq(prod(dims)), ncol = dims[2], byrow = byrow)
    pos <- as.numeric(which(t(lyt) == loc, arr.ind = TRUE))
    plims <- get_plims(dims)
    c(plims$horiz[pos[1]], plims$horiz[pos[1] + 1], plims$vert[pos[2] + 1], plims$vert[pos[2]])
}
get_inset <- function(figlims, subdims) {
    wd <- diff(figlims[1:2])
    ht <- diff(figlims[3:4])
    figlims[c(1, 1, 3, 3)] + subdims*c(wd, wd, ht, ht)
}


dolphin <- readRDS("../data/dolphin.rds")
dolphin.N <- vcount(dolphin)
dolphin.k <- degree(dolphin)
dolphin.kposs <- seq(1, max(dolphin.k))

ba <- readRDS("../data/ba.rds")
ba.N <- vcount(ba)
ba.k <- degree(ba)
ba.kposs <- seq(1, max(ba.k))

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
ht <- 9
wd <- 14
labelsize <- 1.4 # 2
if(save_plots) {
    pdf("../img/degree-sequences.pdf", height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}

dims <- c(5, 2)
subdims <- c(0.4, 0.9, 0.4, 0.9)
mar <- c(2, 2, 1.25, 0) + 0.25
insetmar <- c(2, 2, 0, 0)
insetylim <- c(0, 2)

plot.new()
par(mar = mar)
sp <- 1

                                        # Dolphin
par(fig = get_ploc(sp, dims), new = TRUE)
barplot(table(dolphin.k), names.arg = dolphin.kposs, col = "gray60",
        cex.axis = labelsize, cex.names = labelsize, cex.lab = labelsize,
        xlab = "", ylab = "Frequency")
mtext("(a)", line = 0.4, adj = -0.05, cex = labelsize)
sp <- sp + 1

for(i in 1:3) {
    par(fig = get_ploc(sp, dims), new = TRUE)
    barplot(dolphin.data[[i]], names.arg = dolphin.kposs, col = seq_len(nrow(dolphin.data[[i]])),
            cex.axis = labelsize, cex.names = labelsize, cex.lab = labelsize, ylim = c(0, 100),
            xlab = "", ylab = "Frequency")
    sp <- sp + 1
}
par(fig = get_ploc(sp, dims), new = TRUE)
barplot(dolphin.data[[4]], names.arg = dolphin.kposs, col = seq_len(nrow(dolphin.data[[4]])),
        cex.axis = labelsize, cex.names = labelsize, cex.lab = labelsize, ylim = c(0, 100),
        xlab = "Degree", ylab = "Frequency")
sp <- sp + 1

                                        # BA
baseseq <- 2:20
insetseq <- 21:length(ba.kposs)
insetnames <- ifelse(insetseq %in% seq(25, length(ba.kposs), by = 5), insetseq, NA)

par(fig = get_ploc(sp, dims), new = TRUE)
barplot(sapply(ba.kposs[baseseq], function(k) sum(ba.k == k)), names.arg = ba.kposs[baseseq], col = "gray60",
        cex.axis = labelsize, cex.names = labelsize, cex.lab = labelsize,
        xlab = "", ylab = "")
mtext("(b)", line = 0.4, adj = -0.05, cex = labelsize)
par(fig = get_inset(get_ploc(sp, dims), subdims), new = TRUE, mar = insetmar)
barplot(sapply(ba.kposs[insetseq], function(k) sum(ba.k == k)),
        names.arg = insetnames, col = "gray60",
        cex.axis = 0.75*labelsize, cex.names = 0.75*labelsize, cex.lab = 0.75*labelsize,
        xlab = "", ylab = "", ylim = insetylim, yaxt = "no")
axis(2, at = 0:2, labels = c(0, 1, 2), cex.axis = 0.75*labelsize)
sp <- sp + 1

for(i in 1:3) {
    par(fig = get_ploc(sp, dims), new = TRUE, mar = mar)
    barplot(ba.data[[i]][, baseseq], names.arg = ba.kposs[baseseq], col = seq_len(nrow(ba.data[[i]])),
            cex.axis = labelsize, cex.names = labelsize, cex.lab = labelsize,
            xlab = "", ylab = "")
    par(fig = get_inset(get_ploc(sp, dims), subdims), new = TRUE, mar = insetmar)
    barplot(ba.data[[i]][, insetseq], names.arg = insetnames, col = seq_len(nrow(ba.data[[i]])),
            cex.axis = 0.75*labelsize, cex.names = 0.75*labelsize, cex.lab = 0.75*labelsize,
            xlab = "", ylab = "", ylim = insetylim, yaxt = "no")
    axis(2, at = 0:2, labels = c(0, 1, 2), cex.axis = 0.75*labelsize)
    sp <- sp + 1
}

par(fig = get_ploc(sp, dims), new = TRUE, mar = mar)
barplot(ba.data[[4]][, baseseq], names.arg = ba.kposs[baseseq], col = seq_len(nrow(ba.data[[4]])),
        cex.axis = labelsize, cex.names = labelsize, cex.lab = labelsize,
        xlab = "Degree", ylab = "")
par(fig = get_inset(get_ploc(sp, dims), subdims), new = TRUE, mar = insetmar)
barplot(ba.data[[4]][, insetseq], names.arg = insetnames, col = seq_len(nrow(ba.data[[i]])),
            cex.axis = 0.75*labelsize, cex.names = 0.75*labelsize, cex.lab = 0.75*labelsize,
            xlab = "", ylab = "", ylim = insetylim, yaxt = "no")
axis(2, at = 0:2, labels = c(0, 1, 2), cex.axis = 0.75*labelsize)

if(save_plots) dev.off()

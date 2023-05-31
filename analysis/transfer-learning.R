library(parallel)
ncores <- detectCores() - 1
RNGkind("L'Ecuyer-CMRG")
library(igraph)

args <- commandArgs(trailingOnly = TRUE)
network <- args[1] #"dolphin"
model <- args[2] # "dw"
outfile <- paste0(network, "-", model, ".pdf")

source("../src/functions.R")
load(paste0("../data/", network, ".rda"))
load(paste0("../data/fullstate-", network, ".rda"))
load(paste0("../data/fullstate-alt-", network, ".rda"))

g <- get(network)
N <- vcount(g)
A <- as_adj(g, "both", sparse = FALSE)
k <- degree(g)
n <- floor(log(N))

bp <- switch(
    model,
    dw = seq(0.001, 1, length.out = lout),
    SIS = seq(0.001, 1, length.out = lout),
    genereg = seq(1, 0.001, length.out = lout),
    mutualistic = seq(2, 0.1, length.out = lout),
    wilsoncowan = seq(0.001, 15, length.out = lout)
)
alt1 <- switch(
    model,
    dw = "dw1",
    SIS = "sis1",
    genereg = "gr1",
    mutualistic = "ms1",
    wilsoncowan = "wc1"
)
alt2 <- switch(
    model,
    dw = "dw2",
    SIS = "sis2",
    genereg = "gr2",
    mutualistic = "ms2",
    wilsoncowan = "wc2"
)

GT <- fullstate[[model]] # full ground truth
gt <- rowMeans(fullstate[[model]]) # mean state

soln <- experiment(n, gt, GT, bp)
vs <- soln$vs

ntrials <- 100
errors.ref <- unlist(
        mclapply(seq(ntrials), function(x) experiment(n, gt, GT, bp)$error, mc.cores = ncores)
)
errors.alt1 <- with(list(y = rowMeans(fullstate_alt[[alt1]]), Y = fullstate_alt[[alt1]]), {
    unlist(mclapply(seq(ntrials), function(x) experiment(n, y, Y, bp)$error, mc.cores = ncores))
    ##mean(replicate(25, experiment(n, y, Y, bp)$error))
})
errors.alt2 <- with(list(y = rowMeans(fullstate_alt[[alt2]]), Y = fullstate_alt[[alt2]]), {
    unlist(mclapply(seq(ntrials), function(x) experiment(n, y, Y, bp)$error, mc.cores = ncores))
    ##mean(replicate(25, experiment(n, y, Y, bp)$error))
})

ehat.ref <- obj_fn(vs, gt, GT, bp, optimize_weights = FALSE)
ehat.alt1 <- obj_fn(
    vs, y = rowMeans(fullstate_alt[[alt1]]), Y = fullstate_alt[[alt1]], bp, optimize_weights = FALSE
)
ehat.alt2 <- obj_fn(
    vs, y = rowMeans(fullstate_alt[[alt2]]), Y = fullstate_alt[[alt2]], bp, optimize_weights = FALSE
)

lims <- c(0.025, 0.975)
quantile(errors.ref, probs = lims)
ehat.ref
with(list(e = ehat.ref, x = errors.ref), (e - mean(x))/sd(x))
quantile(errors.alt1, probs = lims)
ehat.alt1
with(list(e = ehat.alt1, x = errors.alt1), (e - mean(x))/sd(x))
quantile(errors.alt2, probs = lims)
ehat.alt2
with(list(e = ehat.alt2, x = errors.alt2), (e - mean(x))/sd(x))

## plotthis <- function(vs, y, Y, bp, legend = FALSE) {
##     palette("Paired")
##     matplot(bp, Y, xlab = "Control parameter", ylab = "x", type = "l",
##             col = 1, lty = 1, lwd = .5,
##             cex.axis = 1.75, cex.lab = 1.75)
##     lines(bp, y, col = 2, lwd = 4, lty = 1)
##     matlines(bp, Y[, vs], lty = 1, lwd = 2, col = 3)
##     lines(bp, rowMeans(Y[, vs]), lty = 1, lwd = 4, col = 4)
##     if(legend) {
##         legend(
##             "topleft", bty = "n", col = 1:4, lty = 1, lwd = 2, cex = 1.5,
##             legend = c("Node states", "Mean state", "Selected nodes", "Approximation")
##         )
##     }
## }

## pdf(outfile, width = 18, height = 6)
## par(mfrow = c(1, 3))
##                                         # original
## plotthis(vs, gt, GT, bp, legend = TRUE)
## ## text(.75*max(bp), .2*max(GT), labels = bquote(r[2] == 3), cex = 1.75, col = "black")
## text(
##     .75*max(bp), .15*max(GT), cex = 1.75, col = "black",
##     labels = bquote(epsilon == .(round(obj_fn(vs, gt, GT, bp, optimize_weights = FALSE), 4)))
## )
##                                         # low
## with(list(y = rowMeans(fullstate_alt[[alt1]]), Y = fullstate_alt[[alt1]]), {
##     plotthis(vs, y, Y, bp)
##     ## text(.75*max(bp), .2*max(Y), labels = bquote(r[2] == 2), cex = 1.75, col = "black")
##     text(
##         .75*max(bp), .15*max(Y), cex = 1.75, col = "black",
##         labels = bquote(epsilon == .(round(obj_fn(vs, y, Y, bp, optimize_weights = FALSE), 4)))
##     )
## })
##                                         # high
## with(list(y = rowMeans(fullstate_alt[[alt2]]), Y = fullstate_alt[[alt2]]), {
##     plotthis(vs, y, Y, bp)
##     ## text(.75*max(bp), .2*max(Y), labels = bquote(r[2] == 4), cex = 1.75, col = "black")
##     text(
##         .75*max(bp), .15*max(Y), cex = 1.75, col = "black",
##         labels = bquote(epsilon == .(round(obj_fn(vs, y, Y, bp, optimize_weights = FALSE), 4)))
##     )
## })
## dev.off()

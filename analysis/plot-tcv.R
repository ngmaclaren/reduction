library(sfsmisc)

network <- "proximity" # "email"
save_plots <- FALSE # TRUE

load(paste0("tcv-", network, ".RData"))

## for (a) and (c)
palette("R4")
ht <- 5
wd <- 5
if(save_plots) {
    pdf(paste0("../img/time-complexity-with-n-", network, ".pdf"), height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
par(mar = c(5, 5, 1, 1), pty = "s")
matplot(
    tc_ns, tc_runtimes, xlab = "", ylab = "", axes = FALSE,
    pch = 1, col = 1, cex = 1.5, ylim = c(0, 1.1*max(tc_runtimes)), xlim = c(0, 31), yaxs = "i", xaxs = "i"
)
box()
eaxis(1, cex.axis = 1.5)
eaxis(2, cex.axis = 1.5)
title(xlab = "n", font.lab = 3, cex.lab = 2)
title(ylab = "Wall time", cex.lab = 2)
if(save_plots) dev.off()


## for (b) and (d)
## str of sapply(results[[1]], get_error): matrix with ns in cols and ntrials in rows
ylim <- range(unlist(lapply(hm_results, function(x) as.numeric(sapply(x, optNS::get_error)))))

## errors <- get_error(results)

palette("Tableau 10")
if(save_plots) {
    pdf(paste0("../img/error-with-hmax-", network, ".pdf"), height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
par(mar = c(5, 5, 1, 1), pty = "s")
plot(NULL, xlim = c(0, 1.05*max(maxits)), ylim = ylim, xlab = "", ylab = "", axes = FALSE, log = "y", xaxs = "i") # xlim = range(maxits)
##     maxits, errors, xlab = "", ylab = "", axes = FALSE, log = "y",
##     cex = 2
## )
for(i in seq_along(maxits)) {
    maxit <- maxits[i]
    for(j in seq_along(hm_ns)) {
        n <- hm_ns[j]
        points(rep(maxit, ntrials), optNS::get_error(hm_results[[i]][[j]]),
               cex = 1, pch = j, col = j)
    }
}
legend("topright", bty = "n", pch = seq_along(hm_ns), col = seq_along(hm_ns),
       legend = hm_ns, title = "n =")
## palette("R4")
abline(v = 50*N, col = "gray50", lwd = 0.5)
box()
eaxis(1, cex.axis = 1.5)
eaxis(2, cex.axis = 1.5, las = 0)
title(xlab = expression(italic(h)[max]), cex.lab = 2)
title(ylab = "Approximation error", cex.lab = 2)
##legend("topright", bty = "n", legend = "50N", lwd = 2, col = 2)
if(save_plots) dev.off()

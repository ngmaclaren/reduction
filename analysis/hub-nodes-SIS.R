library(igraph)
library(sfsmisc)
library(optNS)

net <- "enron"
dyn <- "SIS"
save_plots <- TRUE # FALSE
                                        # network
g <- readRDS(paste0("../data/", net, ".rds"))
                                        # hub nodes
N <- vcount(g)
n <- floor(log(N))
k <- degree(g)
hubs_k <- sort(k, decreasing = TRUE)[1:n]
hubs <- which(k %in% hubs_k)
                                        # node set
ns <- readRDS(paste0("../data/ns-", net, "_", dyn, ".rds"))$opt
best <- ns[[which.min(get_error(ns))]]
                                        # fullstate
Y <- readRDS(paste0("../data/fullstate-", net, ".rds"))[[dyn]]
                                        # approximations
y <- rowMeans(Y)
z <- rowMeans(Y[, best$vs])
zhubs <- rowMeans(Y[, hubs])
                                        # x axis
Ds <- seq(0, 1, length.out = 100)

                                        # plot
if(save_plots) {
    pdf("../img/hub-nodes-SIS.pdf")
} else {
    dev.new()
}
par(mar = c(5, 5, 1, 1), pty = "s")
meancolor <- "black"
sentinelcolor <- "#dc8add"
hubcolor <- "#ed333b"
plot(NULL, xlim = range(Ds), ylim = range(c(y, z, zhubs)), xlab = "D", ylab = "x", font.lab = 3, cex.lab = 2,
     axes = FALSE)
lines(Ds, y, lwd = 6, col = meancolor)
lines(Ds, z, lwd = 3, col = sentinelcolor)
lines(Ds, zhubs, lwd = 3, col = hubcolor)
box()
eaxis(1, cex.axis = 2)
eaxis(2, cex.axis = 2)
legend("bottomright", bty = "n", lwd = 2, col = c(meancolor, sentinelcolor, hubcolor), cex = 1.5,
       legend = c("Average state", "Sentinel node approximation", "Hub node approximation"))
if(save_plots) dev.off()

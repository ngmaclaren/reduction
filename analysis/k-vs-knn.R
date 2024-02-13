## if(interactive()) {
##     if(getwd() != "/user/neilmacl/Documents/reduction/analysis/") {
##         setwd("/user/neilmacl/Documents/reduction/analysis/")
##     }
## }

library(igraph)

network <- "proximity"
dynamics <- "dw"

filetag <- paste(c(network, dynamics), collapse = "-")
infile <- paste0("../data/optimized-nodesets/", filetag, ".rds")
opts <- readRDS(infile)

load(paste0("../data/", network, ".rda"))
g <- get(network)
N <- vcount(g)
n <- floor(log(N))
k <- degree(g)
knn <- knn(g)$knn

vs <- sapply(opts, `[[`, "vs")
ks <- k[vs]
knns <- knn[vs]

pdf(paste0("../img/k-knn-scatter-test-", network, ".pdf"))
plot(k, knn, xlab = "k", ylab = "knn", col = 1, cex = 1, pch = 16)
points(ks, knns, col = adjustcolor(2, .2), pch = 16, cex = 2)
legend("topright", bty = "n", col = 1:2, pch = 16, legend = c("Original network", "Optimized node sets"))
dev.off()

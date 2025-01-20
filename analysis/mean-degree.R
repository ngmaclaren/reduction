library(igraph)
library(sfsmisc)
library(optNS)

dynamics <- "genereg" # doublewell mutualistic SIS genereg

nets <- c(
    "dolphin", "proximity", "celegans", "euroroad", "email",
    "drosophila", "reactome", "route_views", "spanish", "foldoc",
    "tree_of_life", "word_assoc", "enron", "marker_cafe", "prosper",
    "er", "ba", "hk", "gkk", "lfr"
)
nodesets <- lapply(nets, function(network) {
    readRDS(paste0("../data/ns-", network, "_", dynamics, ".rds"))$opt
})
networks <- lapply(nets, function(network) {
    readRDS(paste0("../data/", network, ".rds"))
})
names(networks) <- names(nodesets) <- nets # names(nodefeatures) <- 

p_from_z <- function(z, tail = c("left", "right", "two")) {
    ## https://stats.stackexchange.com/questions/267192/doubling-or-halving-p-values-for-one-vs-two-tailed-tests/267262#267262
    ## Method at that link is more general, but for convenience focusing on z-scores.
    tail <- match.arg(tail)

    fz <- pnorm(z) # this is the default, assuming standard normal (hence z) and integrating from the left
    pos <- z > 0

    switch(
        tail,
        left = fz,
        right = 1 - fz,
        two = switch(
            pos + 1, # make FALSE | TRUE (= 0 | 1) into 1 | 2
            2*fz, # 
            2*(1 - fz)
        )
    )
}

analyze <- function(ns, g, lower.tail = FALSE) { # needs mapply over nodesets and networks
    ks <- get_ks(ns)
    kbar <- mean(degree(g))

    m <- mean(ks)
    
    sd.all <- sd(ks)
    z.all <- (m - kbar)/sd.all
    p.all <- p_from_z(z.all, "two")

    sd.ns <- sd(colMeans(ks))
    z.ns <- (m - kbar)/sd.ns
    p.ns <- p_from_z(z.ns, "two")

    return(data.frame(
        m = m, kbar = kbar,
        std1 = sd.all, z1 = z.all, p1 = p.all,
        std2 = sd.ns, z2 = z.ns, p2 = p.ns
    ))
}

print(dynamics)
(df <- do.call(rbind, mapply(analyze, nodesets, networks, SIMPLIFY = FALSE)))

t.test(df$kbar, df$m, paired = TRUE)

par(mar = c(5, 5, 1, 1))
plot(df$kbar, df$m, xlab = expression(italic(bar(k))), ylab = "m", font.lab = 3, cex.lab = 2, cex = 2,
     axes = FALSE)
abline(a = 0, b = 1, lty = 2, col = 2, lwd = 2)
box()
eaxis(1, cex.axis = 2)
eaxis(2, cex.axis = 2)



### All optimized node sets have the same mean degree for several networks, genereg dynamics
### proximity, foldoc, word_assoc, ER, BA, HK, LFR
### checking a new node set with proximity
network <- "proximity"
dynamics <- "SIS"
fullstatefile <- paste0("../data/fullstate-", network, ".rds")

g <- networks[[network]]
N <- vcount(g)

Y <- readRDS(fullstatefile)[[dynamics]]
y <- rowMeans(Y)
n <- floor(log(N))

soln <- select_optimized(n, g, y, Y)
show(soln)
print(mean(soln$ks)) # 13.5. Same.

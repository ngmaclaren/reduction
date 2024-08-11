library(parallel)
ncores <- detectCores()-1
library(igraph)

save_plots <- TRUE # FALSE

networks <- c( # only the empirical networks
    "dolphin", "celegans", "proximity", "euroroad", "email"
)
graphlist <- lapply(networks, function(network) readRDS(paste0("../data/", network, ".rds")))
dynamics <- c("doublewell", "SIS", "mutualistic", "genereg")
ns.types <- c("opt", "rand")
conds <- expand.grid(networks, dynamics, ns.types, stringsAsFactors = FALSE)
colnames(conds) <- c("networks", "dynamics", "ns.type")
nslistnames <- apply(conds[, 1:2], 1, paste, collapse = "_")
fullstates <- lapply(networks, function(network) readRDS(paste0("../data/fullstate-", network, ".rds")))
names(fullstates) <- networks
nodesets <- lapply(nslistnames, function(cond) {
    readRDS(paste0("../data/ns-", cond, "_w.rds"))
}) 
names(nodesets) <- nslistnames

ks <- lapply(graphlist, degree)
knns <- lapply(graphlist, function(g) knn(g)$knn)
lcls <- lapply(graphlist, transitivity, type = "localundirected")
partitions <- lapply(graphlist, function(g) membership(cluster_louvain(g)))
names(ks) <- names(knns) <- names(lcls) <- names(partitions) <- networks

collector <- function(row) {
    network <- row[1]
    dynamics <- row[2]
    ns.type <- row[3]
    
    ns <- nodesets[[paste(c(network, dynamics), collapse = "_")]][[ns.type]]
    vs <- as.numeric(sapply(ns, `[[`, "vs"))
    ws <- as.numeric(sapply(ns, `[[`, "ws"))
    k <- ks[[network]][vs]
    knn <- knns[[network]][vs]
    lcl <- lcls[[network]][vs]
    lcl[!is.finite(lcl)] <- 0
    comm <- partitions[[network]][vs]
    comm_size <- as.numeric(table(partitions[[network]])[comm])

    data.frame(
        v = vs,
        w = ws, k = k, knn = knn, lcl = lcl, size = comm_size,
        network = network, dynamics = dynamics, ns.type = ns.type,
        row.names = NULL
    )
}

dflist <- apply(conds, 1, collector, simplify = FALSE)
df <- do.call(rbind, dflist)


library(betareg)
## https://stats.stackexchange.com/questions/31300/dealing-with-0-1-values-in-a-beta-regression
df$w. <- (df$w*(nrow(df) - 1) + 0.5)/nrow(df)
df$network <- factor(df$network)
df$dynamics <- factor(df$dynamics)
df$ns.type <- factor(df$ns.type)
df$ns.type <- relevel(df$ns.type, "rand")

model <- betareg(w. ~ poly(k, 2) + knn + lcl + size + network + dynamics + ns.type, data = df)
summary(model)


model.gamma <- glm(
    w. ~ poly(k, 2) + knn + lcl + size + network + dynamics + ns.type,
    data = df,
    family = Gamma
)
summary(model.gamma)





palette("Dark 2")
labelsize <- 2
ht <- 8
wd <- 8
if(save_plots) {
    pdf("../img/ks-ws.pdf", height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
par(mar = c(5, 5, 1, 1), mfrow = c(3, 2))
for(i in seq_along(networks)) {
    net <- networks[i]
    dyn <- "doublewell"
    plotdf <- df[df$network == net & df$dynamics == dyn, ]
    plot(
        x = plotdf$k, xlab = "k",
        y = plotdf$w, ylab = "w",
        type = "p",
        pch = ifelse(as.character(plotdf$ns.type) == "opt", 1, 3),
        col = ifelse(as.character(plotdf$ns.type) == "opt", 1, adjustcolor(3, 0.5)),
        cex.axis = labelsize, cex.lab = labelsize, font.lab = 3
    )
    mtext(paste0("(", letters[i], ")"), adj = 0.05, line = -2.2, cex = labelsize)
    if(i == 1) {
        legend(
            "topright", bty = "n", col = c(3, 1), pch = c(3, 1), cex = 0.75*labelsize, pt.cex = 2, pt.lwd = 2,
            legend = c("Random", "Optimized")
        )
    }
}
if(save_plots) dev.off()

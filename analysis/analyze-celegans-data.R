if(getwd() != "/home/neil/Documents/reduction/analysis/") {
    setwd("/home/neil/Documents/reduction/analysis/")
}

network <- "celegans" # "proximity"

library(igraph)
load(paste0(network, "-data.rda"))
g <- get(network)

N <- vcount(g)
k <- degree(g)
knn <- knn(g)$knn
lcl <- transitivity(g, "localundirected")
bc <- betweenness(g, directed = FALSE)

palette("Tableau 10")
dev.new(width = 14)
plot(NULL, xlim = c(1, N), ylim = c(0.001, 0.7), xlab = "Node (index)", ylab = "p", axes = FALSE, log = "y")
axis(1)
axis(2)
box()
legend("topright", bty = "n", col = 1:3, pch = 16, legend = c("dw", "mut", "sis"))
datas <- list(dw.data, mut.data, sis.data)
for(i in seq_along(datas)) {
    idx <- which(datas[[i]]$p > 0)
    points(c(1:N)[idx], datas[[i]]$p[idx], col = i, pch = 1)
}

## x axis is knn, lcl, bc
## y axis is p
## color is k
cr <- colorRamp(rev(c("red", "orange", "yellow", "green", "blue", "purple")))(sort(unique(log(k)))/max(log(k)))
colors <- apply(cr, 1, function(row) rgb(row[1], row[2], row[3], maxColorValue = 255))
names(colors) <- sort(unique(k))
dev.new()
plot(knn, sis.data$p, col = colors[paste(k)])
## legend("topleft", col = colors, legend = sort(unique(k)), pch = 1)

plot(data.frame(doublewell = dw.data$p, mutualistic = mut.data$p, SIS = sis.data$p), main = "p")
cor(data.frame(dw.data$p, mut.data$p, sis.data$p), method = "spearman")

plot(data.frame(doublewell = dw.data$e, mutualistic = mut.data$e, SIS = sis.data$e), main = "e")
cor(data.frame(dw.data$e, mut.data$e, sis.data$e), method = "spearman")

plot(dw.data$e, sis.data$e)

points(1:N, dw.data$p, col = 1, pch = 1)
points(1:N, mut.data$p, col = 2, pch = 1)
points(1:N, sis.data$p, col = 3, pch = 1)


plot(data.frame(dw.data$e, k, bc, knn, lcl), log = "xy")


plot(data.frame(B.e, k, bc, knn, lcl), log = "xy")

plot(A.e, B.e, log = "xy")


idx <- which(k > 2 & k < 11)
plot(data.frame(A.e, k, bc, knn, lcl)[idx, ], log = "xy")
plot(data.frame(B.e, k, bc, knn, lcl)[idx, ], log = "xy")

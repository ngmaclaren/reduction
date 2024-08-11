library(igraph)
library(optNS)

network <- "lfr" # ba

g <- readRDS(paste0("../data/", network, ".rds"))
N <- vcount(g)
k <- degree(g)
knn <- knn(g)$knn
knnk <- knn(g)$knnk

dynamics <- c("doublewell", "SIS", "mutualistic", "genereg")
names.ns.all <- apply(expand.grid(network, dynamics), 1, function(row) paste(row, collapse = "_"))
nodesets <- lapply(names.ns.all, function(nsname) readRDS(paste0("../data/ns-", nsname, ".rds"))[["opt"]])

allvs <- sapply(nodesets, get_vs)
freqs <- table(as.numeric(allvs))

dev.new(width = 14)
barplot(freqs, names.arg = names(freqs))
chisq.test(freqs)

dev.new()
plot(k[as.numeric(names(freqs))], as.numeric(freqs), xlab = expression(k[i]), ylab = "Frequency", type = "p")

select <- as.numeric(names(freqs)[which.max(freqs)])

k[select]
summary(k)

knn[select]
summary(knn[k == k[select]])

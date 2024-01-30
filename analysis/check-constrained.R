if(interactive()) {
    if(getwd() != "/user/neilmacl/Documents/reduction/analysis/") {
        setwd("/user/neilmacl/Documents/reduction/analysis/")
    }
}

library(igraph)

## didn't "work":
## mean was higher but failed t-test:
## proximity/mutualistic, proximity/genereg
## mean not higher, or really close:
## er/(dw mutualistic genereg); euroroad/dw; email/(dw genereg)
## so there are many, even at 95%. 
network <- "proximity"
dynamics <- "SIS"

filetag <- paste(c(network, dynamics), collapse = "-")
infile <- paste0("../data/optimized-nodesets/", filetag, ".rds")
opts <- readRDS(infile)

load(paste0("../data/", network, ".rda"))
g <- get(network)
N <- vcount(g)
n <- floor(log(N))
k <- degree(g)
knn <- knn(g)$knn

vs <- sapply(opts, `[[`, "vs") # this is a matrix already
ks <- k[vs]
knns <- knn[vs]

top5k <- which(k > quantile(k, .95))

top5 <- as.numeric(V(g)[top5k])
bottom95 <- as.numeric(V(g)[-top5k])

top5counts <- sapply(top5, function(i) sum(apply(vs, 2, function(col) i %in% col)))/ncol(vs)
bottom95counts <- sapply(bottom95, function(i) sum(apply(vs, 2, function(col) i %in% col)))/ncol(vs)

summary(top5counts)
summary(bottom95counts)
t.test(top5counts, bottom95counts, "less")


remaining <- as.numeric(V(g)[-top5k])
rk <- k[remaining]
rknn <- knn[remaining]

top5knn <- which(rknn > quantile(rknn, .95))
bottom5knn <- which(rknn < quantile(rknn, .05))

top5 <- as.numeric(V(g)[top5knn])
bottom5 <- as.numeric(V(g)[bottom5knn])
middle90 <- as.numeric(V(g)[-c(top5knn, bottom5knn)])

top5counts <- sapply(top5, function(i) sum(apply(vs, 2, function(col) i %in% col)))/ncol(vs)
bottom5counts <- sapply(bottom5, function(i) sum(apply(vs, 2, function(col) i %in% col)))/ncol(vs)
middle90counts <- sapply(middle90, function(i) sum(apply(vs, 2, function(col) i %in% col)))/ncol(vs)

summary(top5counts)
summary(bottom5counts)
summary(middle90counts)

tdf <- data.frame(
    counts = c(top5counts, bottom5counts, middle90counts),
    type = c(
        rep("top5", length(top5counts)),
        rep("bottom5", length(bottom5counts)),
        rep("middle90", length(middle90counts))
    )
)
tdf$type <- factor(tdf$type)
tdf$type <- relevel(tdf$type, "middle90")

m <- aov(counts ~ type, data = tdf)
anova(m)
TukeyHSD(m)

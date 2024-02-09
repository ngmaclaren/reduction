## if(interactive()) {
##     if(getwd() != "/user/neilmacl/Documents/reduction/analysis/") {
##         setwd("/user/neilmacl/Documents/reduction/analysis/")
##     }
## }

library(optparse)
optionlist <- list(
    make_option(
        "--network", type = "character", default = "dolphin",
        help = "Default is %default. Options are: dolphin, celegans, proximity, euroroad, email, er, ba, hk, gkk, lfr."
    ),
    make_option(
        "--dynamics", type = "character", default = "dw",
        help = "Default is %default. Options are: dw, SIS, mutualistic, genereg."
    )
)
args <- parse_args(
    OptionParser(option_list = optionlist),
    convert_hyphens_to_underscores = TRUE
)

library(igraph, warn.conflicts = FALSE)

## network <- "celegans"
## dynamics <- "dw"
network <- args$network
dynamics <- args$dynamics
outfile <- "check-constrained.out"

sink(outfile, TRUE, type = "output")
cat("\n", network, dynamics, "\n\n")
sink()

filetag <- paste(c(network, dynamics), collapse = "-")
infile <- paste0("../data/optimized-nodesets/", filetag, ".rds")
opts <- readRDS(infile)

load(paste0("../data/", network, ".rda"))
g <- upgrade_graph(get(network))
N <- vcount(g)
n <- floor(log(N))
k <- degree(g)
knn <- knn(g)$knn

vs <- sapply(opts, `[[`, "vs") # this is a matrix already
ks <- k[vs]
knns <- knn[vs]

top5k <- which(k > quantile(k, 0.90))#.95))

top5 <- as.numeric(V(g)[top5k])
bottom95 <- as.numeric(V(g)[-top5k])

top5counts <- sapply(top5, function(i) sum(apply(vs, 2, function(col) i %in% col)))/ncol(vs)
bottom95counts <- sapply(bottom95, function(i) sum(apply(vs, 2, function(col) i %in% col)))/ncol(vs)

## cat(
##     show(summary(top5counts)),
##     show(summary(bottom95counts))
## )

sink(file = outfile, TRUE, "output")
show(t.test(top5counts, bottom95counts, "less"))
sink()


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

## cat(
##     show(summary(top5counts)),
##     show(summary(bottom5counts)),
##     show(summary(middle90counts))
## )

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

sink(outfile, TRUE, "output")
show(anova(m))
show(TukeyHSD(m))
sink()

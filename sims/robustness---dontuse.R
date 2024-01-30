library(parallel)
ncores <- detectCores()-1
library(igraph)

load("../data/celegans.rda")
load("../data/fullstate.rda")

g <- celegans
N <- vcount(g)
k <- degree(g)
knn <- knn(g)$knn
n <- floor(log(N))

ntrials <- ncores
bparam <- switch()

ec <- eigenvector_centrality(g)
ec_demo ## grab the top n from ec
dc_demo ## grab the top n from k --- are these the same in undirected networks?

dw.env <- new.env()
SIS.env <- new.env()

with(dw.env, {
    dynamics <- "dw"
    bparam <- bparam[]
    Y <- full
    y <- rowMeans(Y)
    
    opt <- make_dataset(n, ntrials, bparam, y, Y, optimize = TRUE)
    bestopt <- opt[which.min(sapply(opt, `[[`, "error"))][[1]]
    fixed <- make_dataset(n, ntrials, bparam, y, Y, comps = bestopt$vs)
    rand <- make_dataset(n, ntrials, bparam, y, Y)
    constr <- make_dataset(n, ntrials, bparam, y, Y, use_connections = TRUE)
    quant <- make_dataset(n, ntrials, bparam, y, Y, use_connections = TRUE, use_quantiles = TRUE)

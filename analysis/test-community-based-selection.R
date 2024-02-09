library(igraph)
source("../src/functions.R")

network <- "celegans"
load(paste0("../data/", network, ".rda"))
g <- upgrade_graph(get(network))
k <- degree(g)
N <- vcount(g)
n <- floor(log(N))
part <- cluster_fast_greedy(g)
sizes <- as.numeric(table(membership(part)))
pvec <- sqrt(sizes)/sum(sqrt(sizes))

select_by_comm_prob <- function(n, g, partition, pvec, vs = V(g)) {
    stopifnot(is.numeric(n))
    stopifnot(is.igraph(g))
    stopifnot(inherits(part, "communities"))
    stopifnot(length(pvec) == length(part))

    withheld <- V(g)[-which(V(g) %in% vs)]
    if(length(withheld) > 0) {
        available <- V(g)[-withheld]
    } else {
        available <- vs
    }
    Cs <- seq_len(length(partition))
    mbr <- membership(partition)
    
    C_vec <- sample(Cs, n, TRUE, pvec) # this tells me which communities to take from

    df <- data.frame(v = as.numeric(available), C = as.numeric(mbr[available]))

    sapply(C_vec, function(C) sample(df$v[df$C == C], 1))
    
}

select_by_community <- function(n, vs, part = NULL, alg = "cluster_louvain") {
    stopifnot(is.numeric(n))
    stopifnot(is.igraph(g))

    if(is.null(part)) {
        part <- get(alg)(g)
    } else {
        stopifnot(inherits(part, "communities"))
    }
    
    stopifnot(length(part) >= n)

    available <- vs
    mbr <- membership(part)
    S <- c()
    m <- c()

    while(length(S) < n) {
        
        v <- sample(available, 1)
        S <- c(S, v)
        m <- c(m, mbr[v])
        
        toremove <- which(mbr %in% m)
        available <- vs[-toremove]
    }

    return(S)
}

## S <- select(n, g, part)
## membership(part)[S]

S <- select_by_community(n, V(g))
membership(part)[S]

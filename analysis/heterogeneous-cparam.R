library(parallel)
ncores <- detectCores()-1
library(igraph)
library(deSolve)
library(sdn)
library(optNS)

network <- "proximity"

dwij <- function(t, x, params) {
    with(params, {
        coupling <- mapply(function(vs, ds) sum(((D*ds)+baseD)*x[vs]), AL, DL) # where D is as usual, AL is the adj list, and DL is a list of indicators \in {0, 1} showing which edges increase in D
        dx <- -(x - r[1])*(x - r[2])*(x - r[3]) + coupling + u
        return(list(c(dx)))
    })
}

genY <- function(condition, prob) { # prob is the proportion of affected (= fixed, don't change)
                                        # Conditions
    E(g)$marked <- switch(
            condition, 
            random = sample(c(0, 1), length(E(g)), replace = TRUE, prob = c(prob, 1 - prob)),
            degsum.high = ifelse(degsum > quantile(degsum, 1 - prob), 0, 1), 
            degsum.low = ifelse(degsum < quantile(degsum, prob), 0, 1),
            edgebet.high = ifelse(edgebet > quantile(edgebet, 1 - prob), 0, 1),
            edgebet.low = ifelse(edgebet < quantile(edgebet, prob), 0, 1)
        )
                                        # For all conditions, make the list
    DL <- lapply(AEL, function(es) E(g)$marked[es])
    params <- c(.doublewell, list(AL = AL, AEL = AEL, DL = DL), baseD = 0.05)
    solve_in_range(Ds, "D", dwij, rep(params$xinit.low, N), params, control, "ode")
}

ns <- readRDS(paste0("../data/ns-", network, "_doublewell.rds")) # hard-coding doublewell
g <- readRDS(paste0("../data/", network, ".rds")) ## on the cluster do proximity and email
N <- vcount(g)

AL <- as_adj_list(g, "all") # returns a list of igraph.vs
AEL <- as_adj_edge_list(g, "all") # returns a list of igraph.es
Ds <- seq(0, 1, length.out = 100)

                                        # degree sum
degsum <- sapply(E(g), function(e) apply(ends(g, e), 1, function(row) sum(degree(g, v = row))))
                                        # edge betweenness
edgebet <- edge_betweenness(g, directed = FALSE)

times <- 0:15
control <- list(times = times, ncores = ncores)

conditions <- c("random", "degsum.high", "degsum.low", "edgebet.high", "edgebet.low")
probs <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
Ylist <- lapply(probs, function(prob) lapply(conditions, genY, prob = prob))
names(Ylist) <- as.character(probs)
for(i in seq_along(Ylist)) names(Ylist[[i]]) <- conditions

## from here, move to a separate plotting file.
rand_error <- get_error(ns$rand)
opt_error <- get_error(ns$opt)
het_errors <- lapply(
    Ylist, function(Y)
        lapply(
            Y, function(Y.)
                sapply(ns$opt, function(x) obj_fn(x$vs, rowMeans(Y.), Y.))
        )
)
xlim <- range(unlist(c(het_errors, rand_error, opt_error)))

save.image(paste0("heterogeneous-cparam-", network, ".RData"))

library(parallel)
ncores <- detectCores()-1
library(igraph)
library(sdn)
library(optNS)

network <- "email"
dynamics <- "doublewell"

fullstatefile <- paste0("../data/fullstate-", network, ".rds")
nsfile <- paste0("../data/ns-", network, "_", dynamics, ".rds")

                                        # this is the original, against which we'll test
g <- readRDS(paste0("../data/", network, ".rds")) 
N <- vcount(g)
M <- ecount(g)

                                        # Make node names to keep track of nodes after subgraph()
V(g)$names <- as.character(as.numeric(V(g)))

                                        # the test dynamics results
Y <- readRDS(fullstatefile)[[dynamics]] 
y <- rowMeans(Y)
Ds <- seq(0, 1, length.out = 100)
n <- floor(log(N))
nodesets <- readRDS(nsfile)[c("opt", "rand")]

control <- list(times = 0:15, ncores = ncores)
ntrials <- 5

### EDGES

## First, remove 10%, 20%, and 30% of edges
edges_toremove <- floor(M * c(0.01, 0.03, seq(0.05, 0.5, by = 0.05))) # remove this many edges from g
removethese_edges <- lapply(edges_toremove, function(x) sample(seq(M), size = x)) # a list of integer vectors. These integer vectors identify edges to remove from g

glistE <- lapply(removethese_edges, function(x) delete_edges(g, x))

                                        # check
## sapply(glist, vcount)
## sapply(glist, ecount)
## dev.new(height = 6, width = 15)
## par(mfrow = c(2, 5), pty = "s")
## for(g. in glist) plot(g., vertex.size = 5, vertex.label = "")

## Now, need ground truth sims for glist
YlistE <- lapply(
    glistE, function(g.) {
        with(list(
            params = c(.doublewell, list(AL = as_adj_list(g., "all")))
        ), {
            solve_in_range(
                params$Ds, "D", doublewell, rep(params$xinit.low, N), params, control, "ode",
                method = "adams", maxsteps = 20000
            )
        })
    }
)

## dev.new(height = 6, width = 15)
## par(mfrow = c(2, 5), pty = "s")
## for(df in Ylist) bifplot(df, Ds, TRUE, lwd = 0.5, col = 1)

uncertainE <- mapply(
    ## function(g., Y.) select_optimized(n, g., rowMeans(Y.), Y.),
    function(g., Y.) {
        make_dataset(
            ntrials = ntrials, ns.type = "opt", ncores = ncores,
            n = n, g = g., y = rowMeans(Y.), Y = Y.
        )
    }, glistE, YlistE, SIMPLIFY = FALSE
) 

uncertain_errorE <- t(sapply(uncertainE, function(unc) sapply(unc, function(ns) obj_fn(ns$vs, y, Y))))


### NODES

## Removing nodes shouldn't be too much harder.
nodes_tokeep <- ceiling(N * c(0.99, 0.97, seq(0.95, 0.5, by = -0.05)))
keepthese_nodes <- lapply(nodes_tokeep, function(x) sample(seq(N), size = x))
glistN <- lapply(keepthese_nodes, function(x) subgraph(g, x))
nsN <- sapply(glistN, function(g.) floor(log(vcount(g.))))

                                        # check
## sapply(glist, vcount)
## sapply(glist, ecount)
## dev.new(height = 6, width = 15)
## par(mfrow = c(2, 5), pty = "s")
## for(g. in glist) plot(g., vertex.size = 5, vertex.label = "")

YlistN <- lapply(
    glistN, function(g.) {
        with(list(
            params = c(.doublewell, list(AL = as_adj_list(g., "all")))
        ), {
            solve_in_range(
                params$Ds, "D", doublewell, rep(params$xinit.low, N), params, control, "ode",
                method = "adams", maxsteps = 20000)
        })
    }
)

## dev.new(height = 6, width = 15)
## par(mfrow = c(2, 5), pty = "s")
## for(df in Ylist) bifplot(df, Ds, TRUE, lwd = 0.5, col = 1)

uncertainN <- mapply(
    ## function(n., g., Y.) select_optimized(n., g., rowMeans(Y.), Y.),
    function(n., g., Y.) {
        make_dataset(
            ntrials = ntrials, ns.type = "opt", ncores = ncores,
            n = n., g = g., y = rowMeans(Y.), Y = Y.
        )
    }, nsN, glistN, YlistN,
    SIMPLIFY = FALSE
)

uncertain_errorN <- t(mapply(
    function(unc, g.) {
        sapply(unc, function(ns) {
            obj_fn(as.numeric(V(g)[which(V(g)$names %in% V(g.)$names[ns$vs])]), y = y, Y = Y)
        })
    },
    uncertainN, glistN
))



save.image(paste0("nu-", network, ".RData"))

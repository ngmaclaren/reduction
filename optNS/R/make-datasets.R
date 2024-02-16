make_dataset <- function(ntrials,
                         ns.type = c("opt", "rand", "fixed", "constr", "quant", "comm"),
                         ncores = 1,
                         ...) {
    ns.type <- match.arg(ns.type) # don't remember how to do this

    selector <- switch(
        ns.type,
        opt = select_optimized,
        rand = select_randomized,
        fixed = select_bydegseq,
        constr = select_constrained,
        quant = select_quantiled,
        comm = select_bycommunity
    )

    with(list(args = list(...)), {
        mclapply(seq_len(ntrials), function(x) do.call(selector, args), mc.cores = ncores)
    })
}


## make_dataset <- function(n, # how many nodes are in the node set
##                          ntrials, # how many times should we do the experiment
##                          bparam, # the bifurcation parameter values (a sequence)
##                          y, # make sure this function accepts the right mean state
##                          Y, # and full state
##                          comps = NULL, # optionally, a degree distribution to fix, provided as vs
##                          optimize = FALSE, # should de novo optimized solutions be produced?
##                          optimize_weights = FALSE, # should weights be optimized?
##                          use_connections = FALSE, # should the ks/knns be used to select nodes?
##                          use_quantiles = FALSE, # if use_connections, also split into n by quantile?
##                          use_communities = FALSE,
##                          verbose = FALSE # for `trace` and `mc.silent`; doesn't seem to do anything
##                          ) {
##     if(is.numeric(comps)) { # If we have some reference node set to match (e.g. a good node set)
##         ks <- k[comps]
##         poss <- lapply(ks, function(ki) as.numeric(V(g)[which(k == ki)]))
##         VS <- replicate(ntrials, sapply(poss, sample, 1), simplify = FALSE) # fixed-degree
##     } else {
##         if(optimize) { # fully optimize
##             dl <- mclapply(
##                 1:ntrials, 
##                 function(x) experiment(
##                                 n, y, Y, bparam, optimize_weights = optimize_weights,
##                                 trace = isTRUE(verbose)
##                             ),
##                 mc.cores = ncores, mc.silent = isFALSE(verbose)#FALSE
##             )
##             return(dl) # can exit with this list of optimized node sets
##         } else if(use_connections) { # impose degree/knn based constraints
##                                         # NB: sample(..., replace=FALSE) is default
##             library(igraph)
##             k <- degree(g)
##             knn <- knn(g)$knn

##                                         # Effect of knn does not appear to be robust
##             ## top10k <- which(k > quantile(k, probs = 0.9))
##             ## bottom10knn <- which(knn < quantile(knn, probs = 0.1))
##             ## avail <- as.numeric(V(g)[-unique(c(top10k, bottom10knn))])
##                                         # Better results with only eliminating the very largest nodes
##             top5k <- which(k > quantile(k, probs = 0.95))
##             avail <- as.numeric(V(g)[-top5k])
            
##             if(use_quantiles) { # impose stratified sampling by degree
##                 probwidth <- 1/n
##                 probs <- seq(from = 0, to = 1, by = probwidth)
##                 quantiles <- quantile(k[avail], probs)
##                 bins <- mapply(
##                     function(x, y) {
##                         if(x == y) {
##                             intersect(which(k == x), avail)
##                         } else {
##                             intersect(which(k >= x & k < y), avail)
##                         }
##                     },
##                     quantiles[-length(quantiles)], quantiles[-1], SIMPLIFY = FALSE
##                 )
                
##                 VS <- replicate(ntrials, sapply(bins, sample, 1), simplify = FALSE) #use_connections & use_quantiles
##             } else if(use_communities) { # impose sampling by community structure
##                 partition <- cluster_louvain(g)
##                                         # use a deterministic algorithm for reproducibility
##                 ## partition <- cluster_fast_greedy(g)

##                 sizes <- as.numeric(table(membership(partition)))
##                 pvec <- sqrt(sizes)/sum(sqrt(sizes))

##                 VS <- replicate(# use_connections & use_communities 
##                     ntrials,
##                     as.numeric(select_by_comm_prob(n, g, partition, pvec, avail)),
##                     simplify = FALSE
##                 )
##             } else {
##                 VS <- replicate(ntrials, sample(avail, n), simplify = FALSE) # use_connections, but no further
##             }
            
##         } else {
##             VS <- replicate(ntrials, sample(1:N, n), simplify = FALSE) # fully random
##         }
##     }

##     KS <- lapply(VS, function(vs) k[vs])
##     sorting <- lapply(KS, order)

##     error <- mclapply(
##         VS, function(x) obj_fn(x, y, Y, bparam, optimize_weights = optimize_weights),
##         mc.cores = ncores, mc.silent = isFALSE(verbose)
##     )

##     dl <- vector("list", length(VS))
##     for(i in seq_along(dl)) {
##         dl[[i]] <- list(
##             vs = VS[[i]][sorting[[i]]],
##             error = error[[i]],
##             ks = KS[[i]][sorting[[i]]]
##         )
##     }
##     if(optimize_weights) {
##         WS <- mclapply(VS, quadoptm, y = y, Y = Y, mc.cores = ncores, mc.silent = isFALSE(verbose))
        
##         for(i in seq_along(dl)) dl[[i]]$ws <- WS[[i]][sorting[[i]]]
##     }

##     return(dl)
## }

## ### Can I get rid of this function?
## generate_data <- function(g, Y, y, n, ntrials = 100, optimize_weights = FALSE, verbose = FALSE) {
##     random <- make_dataset(n, ntrials, doublewell_parms$Ds, y, Y)
##     random_errors <- sapply(random, `[[`, "error")
##     optimized <- make_dataset(
##         n, ntrials, doublewell_parms$Ds, y, Y, optimize = TRUE, optimize_weights = optimize_weights,
##         verbose = verbose
##     )
##     optimized_errors <- sapply(optimized, `[[`, "error")
##     idx <- which.min(optimized_errors)
##     ref <- optimized[[idx]]$vs
##     fixingdeg <- make_dataset(n, ntrials, doublewell_parms$Ds, y, Y, comps = ref)
##     fixingdeg_errors <- sapply(fixingdeg, `[[`, "error")

##     rvs <- t(sapply(random, `[[`, "vs"))
##     rks <- t(sapply(random, `[[`, "ks"))
##     fvs <- t(sapply(fixingdeg, `[[`, "vs"))
##     fks <- t(sapply(fixingdeg, `[[`, "ks"))
##     ovs <- t(sapply(optimized, `[[`, "vs"))
##     oks <- t(sapply(optimized, `[[`, "ks"))

##     dl <- list(
##         re = random_errors,
##         fe = fixingdeg_errors,
##         oe = optimized_errors,
##         rvs = rvs,
##         fvs = fvs,
##         ovs = ovs,
##         rks = rks,
##         fks = fks,
##         oks = oks
##     )
##     return(dl)
## }

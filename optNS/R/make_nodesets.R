## This is a file for the different kinds of random node selection
## They will feed into the make_dataset() function, or be used alone.
## Each needs to return something like this:
## return(list(vs = vs[sorting], error = error, ks = ks[sorting]))

## arguments:
## - n, y, Y
## - ntrials
## - optimize_weights ---> separate this into a different function

## need to pass g in at the top level

## how long does quadoptm() take to do once? Is it a problem to do it again (I mean, after optimizing based on weight, is it the same answer every time to give it the same sequence of nodes? Should be.)

## Problem: some of these functions return lists, others return a single node set
## Solution: all the select_* functions should return a single node set

                                        # Remove Ds arg?
select_optimized <- function(n, g, y, Y, Ds, maxit = NULL, optimize_weights = FALSE,
                             trace = FALSE, sorted = TRUE) {
    if(is.null(maxit)) {
        maxit <- switch(
            optimize_weights+1, # convert from {0, 1} to {1, 2}
            50*ncol(Y), # more iterations if not optimizing weights
            25*ncol(Y)) # fewer iterations if optimizing weights
    }

    k <- degree(g)
    vs <- sample(as.numeric(V(g)), n)
    result <- optim(
        vs, obj_fn, update_vs, y = y, Y = Y, Ds = Ds, # replace with g?
        optimize_weights = optimize_weights,
        method = "SANN", control = list(trace = trace, maxit = maxit, temp = 10)
    )

    vs <- result$par
    
    make_dl(vs, y, Y, k, optimize_weights, sorted)
}

#' Fully random node set construction
#' @param ntrials
#' @param n
#' @param N
#' @return A vector of nodes
select_randomized <- function(n, g, y, Y, optimize_weights = FALSE, sorted = TRUE) {
    k <- degree(g)

    vs <- sample(as.numeric(V(g)), n)

    make_dl(vs, y, Y, k, optimize_weights, sorted)
}

## select_randomly <- function(ntrials, n, N) {
##     replicate(ntrials, sample(1:N, n), simplify = FALSE)
## }

#' Constrained by degree
#' @param k The vector of node degrees
select_constrained <- function(n, g, y, Y, optimize_weights = FALSE, sorted = TRUE) {
    k <- degree(g)
    top5k <- which(k > quantile(k, probs = 0.95))
    avail <- as.numeric(V(g)[-top5k])

    vs <- sample(avail, n)
    
    make_dl(vs, y, Y, k, optimize_weights, sorted)
}

#' Quantiles of the degree distribtuion
#' @param k The vector of node degrees
select_quantiled <- function(n, g, y, Y, optimize_weights = FALSE, sorted = TRUE, constrain = TRUE) {
    k <- degree(g)
    
    if(constrain) {
        top5k <- which(k > quantile(k, probs = 0.95))
        avail <- as.numeric(V(g)[-top5k])
    } else {
        avail <- as.numeric(V(g))
    }

    probwidth <- 1/n
    probs <- seq(from = 0, to = 1, by = probwidth)
    quantiles <- quantile(k[avail], probs)
    bins <- mapply(
        function(x, y) {
            if(x == y) {
                intersect(which(k == x), avail)
            } else {
                intersect(which(k >= x & k < y), avail)
            }
        },
        quantiles[-length(quantiles)], quantiles[-1], SIMPLIFY = FALSE
    )
    
    vs <- sapply(bins, sample, 1)

    make_dl(vs, y, Y, k, optimize_weights, sorted)
}

select_by_comm_prob <- function(n, g, partition, pvec, vs = V(g)) {
                                        # This function selects nodes according to a probability that is
                                        # proportional to the square root of community size.
    stopifnot(is.numeric(n))
    stopifnot(is.igraph(g))
    stopifnot(inherits(partition, "communities"))
    stopifnot(length(pvec) == length(partition))

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

#' Based on community membership
#' @param partition An igraph communities object
#'
select_bycommunity <- function(n, g, y, Y, optimize_weights = FALSE, sorted = TRUE, constrain = TRUE) {
    k <- degree(g)
    
    if(constrain) {
        top5k <- which(k > quantile(k, probs = 0.95))
        avail <- as.numeric(V(g)[-top5k])
    } else {
        avail <- as.numeric(V(g))
    }

    partition <- cluster_louvain(g)
    sizes <- as.numeric(table(membership(partition)))
    pvec <- sqrt(sizes)/sum(sqrt(sizes))

    vs <- as.numeric(select_by_comm_prob(n, g, partition, pvec, avail))

    make_dl(vs, y, Y, k, optimize_weights, sorted)
}

#' Based on a fixed (usually optimal) degree sequence
#' @param comps
select_bydegseq <- function(n, g, comps, y, Y, optimize_weights = FALSE, sorted = TRUE) {
    k <- degree(g)
    ks <- k[comps]
    poss <- lapply(ks, function(ki) as.numeric(V(g)[which(k == ki)]))
    vs <- sapply(poss, sample, 1)

    make_dl(vs, y, Y, k, optimize_weights, sorted)
}


make_dl <- function(vs, y, Y, bparam, k, optimize_weights = FALSE, ws = NULL, sorted = TRUE) {
    if(optimize_weights & is.null(ws)) ws <- quadoptm(vs, y, Y)
    error <- obj_fn(vs, y, Y, bparam, optimize_weights = optimize_weights, ws = ws) 
    ks <- degree(g)[vs]

    if(sorted) {
        sorting <- order(ks)
        vs <- vs[sorting]
        ks <- ks[sorting]
        if(optimize_weights) ws <- ws[sorting]
    }
    
    return(list(vs = vs, error = error, ks = ks, ws = ws))
}
        

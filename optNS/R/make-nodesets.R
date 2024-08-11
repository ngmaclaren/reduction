optim <- stats::optim
quantile <- stats::quantile

sample.local <- function(x, ...) {
    ## sample() can have surprising behavior when sample(x, ...) and x can have varying lengths
    ## this function should remove that surprise
    ## see the example in ?sample
    x[sample.int(length(x), ...)]
}

sample.frombins <- function(bins) {
    ## should return one element from each bin with no repetition
    sample.frombin <- function(bin, out) {
        avail <- bins[[i]][!(bins[[i]] %in% out)]
        if(length(avail) == 0) {
            return(NA)
        } else {
            return(sample.local(avail, size = 1))
        }
    }
    out <- numeric(length(bins))
    for(i in seq_along(bins)) out[i] <- sample.frombin(bins[[i]], out)

    return(out)
}


#' Functions to select sentinel node sets
#'
#' @description Each of these functions returns a node set along with some additional information. See details.
#'
#' @param n The number of nodes to use in the node set.
#' @param g The network from which the nodes will be drawn. 
#' @param y The mean state of the nodes on the network at each value of the bifurcation parameter. 
#' @param Y The state of all nodes at each value of the bifurcation parameter. 
#' @param optimize_weights Whether or not to optimize node weights by quadratic programming.
#' @param sorted Whether or not to sort the nodes by degree in the returned list.
#' @param maxit The default is typically sufficient, but can be set here. 
#' @param trace Whether or not to observe the progress of optim().
#' @param constrain Whether or not to remove the top 5\% largest nodes by degree
#' @param comps The node index sequence
#' @param vs A node set (numeric node indices)
#' @param k The degree distribution of a network; a numeric vector
#' @param ws A vector of node weights
#'
#' @name selector
#'
#' @details These functions all accept similar arguments and return the same data structure (i.e., n node indices, the associated approximation error, etc) but differ in the algorithm used to select the nodes. No algorithm accepts repeated nodes in the output. 
#' - select_optimized() uses combinatorial simulated annealing
#' - select_randomized() selects nodes uniformly at random
#' - select_constrained() removes the top 5\% of nodes by degree from consideration, then selects randomly.
#' - select_knn_constrained() removes the bottom 5\% by average nearest neighbor degree and optionally the top 5\% by degree, then selects randomly
#' - select_quantiled() divides the degree distribution into n bins and selects one node at random from each bin
#' - select_bycommunity() assigns a weight to each node to encourage, but not require, selecting nodes from different communities.
#' - select_bydegseq() accepts a set of nodes as an argument (`comps`) and returns a set of nodes with the same degree sequence, but which is otherwise random.
#'
#' make_dl() accepts a node set `vs` and returns the standard data structure.
#'
#' @return A list with the following elements: vs: The node set, a numeric vector of node indices; error: The approximation error for the node set; ks: The degree sequence; ws: NULL, or the weights assigned to each node if requested.
NULL

#' @rdname selector
#' @export
select_optimized <- function(n, g, y, Y,
                             optimize_weights = FALSE, sorted = TRUE, 
                             maxit = NULL, trace = FALSE) {
    if(is.null(maxit)) {
        maxit <- switch(
            optimize_weights+1, # convert from {0, 1} to {1, 2}
            50*ncol(Y), # more iterations if not optimizing weights
            25*ncol(Y)) # fewer iterations if optimizing weights
    }

    k <- degree(g)
    vs <- sample.local(as.numeric(V(g)), n)
    result <- optim(
        par = vs, fn = obj_fn, gr = update_vs, g = g, y = y, Y = Y,
        optimize_weights = optimize_weights,
        method = "SANN", control = list(trace = trace, maxit = maxit, temp = 10)
    )

    vs <- result$par
    
    make_dl(vs, g, y, Y, k, optimize_weights, sorted = sorted)
}

#' @rdname selector
#' @export
select_randomized <- function(n, g, y, Y, optimize_weights = FALSE, sorted = TRUE) {
    k <- degree(g)

    vs <- sample.local(as.numeric(V(g)), n)

    make_dl(vs, g, y, Y, k, optimize_weights, sorted = sorted)
}

#' @rdname selector
#' @export
select_constrained <- function(n, g, y, Y, optimize_weights = FALSE, sorted = TRUE) {
    k <- degree(g)
    top5k <- which(k > quantile(k, probs = 0.95))
    avail <- as.numeric(V(g)[-top5k])

    vs <- sample.local(avail, n)
    
    make_dl(vs, g, y, Y, k, optimize_weights, sorted = sorted)
}

#' @rdname selector
#' @export
select_knn_constrained <- function(n, g, y, Y, optimize_weights = FALSE, sorted = TRUE, constrain = TRUE) {
    k <- degree(g)
    knn <- knn(g)$knn
    bottom5knn <- which(knn < quantile(knn, probs = 0.05))
    
    if(constrain) {
        top5k <- which(k > quantile(k, probs = 0.95))
        avail <- as.numeric(V(g)[-unique(c(top5k, bottom5knn))])
    } else {
        avail <- as.numeric(V(g)[-bottom5knn])
    }

    vs <- sample.local(avail, n)

    make_dl(vs, g, y, Y, k, optimize_weights, sorted = sorted)
}
    

#' @rdname selector
#' @export
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
    quantiles <- quantile(k[avail], probs, names = FALSE)
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
    
    ## vs <- sapply(bins, sample.local, 1)
    vs <- sample.frombins(bins)

    make_dl(vs, g, y, Y, k, optimize_weights, sorted = sorted)
}

select_by_comm_prob <- function(n, g, partition, pvec, vs = V(g)) {
                                        # This function selects nodes according to a probability that is
                                        # proportional to the square root of community size.
    stopifnot(is.numeric(n))
    stopifnot(is_igraph(g))
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
    
    C_vec <- sample.local(Cs, size = n, replace = TRUE, prob = pvec) # this tells me which communities to take from

    df <- data.frame(v = as.numeric(available), C = as.numeric(mbr[available]))

    ## sapply(C_vec, function(C) sample.local(df$v[df$C == C], 1))
    bins <- lapply(C_vec, function(C) df$v[df$C == C])
    sample.frombins(bins)
    
}

#' @rdname selector
#' @export
select_bycommunity <- function(n, g, y, Y, optimize_weights = FALSE, sorted = TRUE, constrain = TRUE) {
    k <- degree(g)
    
    if(constrain) {
        top5k <- which(k > quantile(k, probs = 0.95))
        avail <- as.numeric(V(g)[-top5k])
    } else {
        avail <- as.numeric(V(g))
    }

    if(is_directed(g)) {
        partition <- cluster_walktrap(g)
    } else {
        partition <- cluster_louvain(g)
    }
    sizes <- as.numeric(table(membership(partition)))
    pvec <- sqrt(sizes)/sum(sqrt(sizes))

    vs <- as.numeric(select_by_comm_prob(n, g, partition, pvec, avail))

    make_dl(vs, g, y, Y, k, optimize_weights, sorted = sorted)
}

#' @rdname selector
#' @export
select_bydegseq <- function(n, g, comps, y, Y, optimize_weights = FALSE, sorted = TRUE) {
    k <- degree(g)
    ks <- k[comps]
    poss <- lapply(ks, function(ki) as.numeric(V(g)[which(k == ki)]))
    ## vs <- sapply(poss, function(vec) sample.local(vec, size = 1)) # allows repetition
    vs <- sample.frombins(poss)

    make_dl(vs, g, y, Y, k, optimize_weights, sorted = sorted)
}

#' @rdname selector
#' @export
make_dl <- function(vs, g, y, Y, k, optimize_weights = FALSE, ws = NULL, sorted = TRUE) {
    if(optimize_weights & is.null(ws)) ws <- quadoptm(vs, y, Y)
    error <- obj_fn(vs, y, Y, optimize_weights = optimize_weights, ws = ws) 
    ks <- degree(g)[vs]

    if(sorted) {
        sorting <- order(ks)
        vs <- vs[sorting]
        ks <- ks[sorting]
        if(optimize_weights) ws <- ws[sorting]
    }
    
    return(list(vs = vs, error = error, ks = ks, ws = ws))
}
        

#' @import ROI ROI.plugin.qpoases
library(ROI)
library(ROI.plugin.qpoases)

### Objective function
## z is the reduction, y is the mean system state
### ^-- both of these are over all values of the bifurcation parameter (i.e., index by l)
## Y is the full matrix of x_{i,l}^*
## This is unfortunately confusing notation, but is used consistently throughout this analysis.
#' Objective function
#'
#' Compute the objective function for optimization.
#'
#' @param z A vector of length L
#' @param y A vector of length L
#'
#' @details
#' This function expects that some mean state, y, can be or has been computed for a system. The function then tests the approximation, z, against the true mean state, y. The objective function is the sum of squared errors between z and y, normalized by the length L times the grand mean (the mean of y). The value that it returns is the approximation error for the given approximation z to the true mean state y. 
#'
#' @return A scalar
#' @examples
#' y <- 1:10
#' z <- sample(1:10, 10)
#' calc_obj(z, y)
#' @export
calc_obj <- function(z, y) sum((z - y)^2)/(length(y)*mean(y))

## #' Compute the reduction of a sample
## #'
## #' Given a set of nodes and the ground truth simulation results, compute either the weighted or unweighted mean state of the sample.
## #' @param vs A sample of nodes (numeric or igraph.vs)
## #' @param y The ground truth mean state
## #' @param Y A matrix
## #' @param Ds The bifurcation parameter. Present for passing arguments inside an optim() call. Can be NULL if the function is used on its own.
## #' @param optimize_weights Logical for whether or not to optimize node weights.
## #' @return The vector of approximated mean states, called `z` in other functions in this package.
## #' 
## #' @export
## reduction.sample <- function(vs, y, Y, Ds,
##                              optimize_weights = FALSE) {
##     Z <- as.matrix(Y[, vs])
##     if(optimize_weights) {
##         w <- quadoptm(vs, y, Y)
##         apply(Z, 1, weighted.mean, w)
##     } else {
##         ##apply(Z, 1, system_state, angle)
##         rowMeans(Z)
##     }
## }

#' Calculate the approximated mean state for a sample of nodes
#'
#' Given a sample of nodes, compute the objective function: the weighted or unweighted mean state of the sample of nodes, normalized appropriately.
#' @param vs A sample of nodes (numeric or igraph.vs)
#' @param y The ground truth mean state
#' @param Y A matrix
#' @param Ds The bifurcation parameter. Present for passing arguments inside an optim() call. Can be NULL if the function is used on its own.
#' @param optimize_weights Logical for whether or not to optimize node weights.
#' 
#' @return A numeric vector of approximated system states at each value of the bifurcation parameter
#' 
#' @seealso [calc_obj()], [reduction.sample()]
#' @export
                                        # can I eliminate `Ds` from the formals? Maybe with `...`?
obj_fn <- function(vs, y, Y, Ds, optimize_weights = FALSE, ws = NULL) { 
    Z <- as.matrix(Y[, vs])

    if(optimize_weights) {
        if(is.null(ws)) ws <- quadoptm(vs, y, Y)
        z <- apply(Z, 1, weighted.mean, ws)
    } else {
        z <- rowMeans(Z)
    }

    calc_obj(z, y)
}


### Optimization
## take_sample <- function(vs, Y) Y[, vs]

#' Quadratic programming function
#'
#' Optimization of node weights by quadratic programming, using ROI and the qpoases algorithm.
#'
#' @param vs The selected nodes (i.e., the columns of Y)
#' @param y The mean state at each level of the bifurcation parameter
#' @param Y The state of each node (in columns) at each level of the bifurcation parameter (rows)
#'
#' @return A vector of weights for the nodes in `vs`.
#' @export
quadoptm <- function(vs, y, Y) {
    n <- length(vs)
    
    Yp <- Y[, vs]
    Q0 <- 2*(t(Yp)%*%Yp) # I had to include the 2, maybe to match the 2 in a0
    a0 <- -2*(t(y)%*%Yp)

    qp <- OP(
        Q_objective(Q0, a0),
        L_constraint(rep(1, n), "==", 1)
    )
    qpsol <- ROI_solve(qp, "qpoases")

    return(solution(qpsol))
}

#' Select a new node
#'
#' Given a starting sample of nodes, select a random node from the node set and replace it with a random node not in the node set.
#' @param vs The selected nodes (i.e., the columns of Y)
#' @param y The mean state at each level of the bifurcation parameter
#' @param Y The state of each node (in columns) at each level of the bifurcation parameter (rows)
#' @param optimize_weights Logical for whether or not to optimize node weights.
#'
#' @return A numeric vector of node indices
#' @export
update_vs <- function(vs, y, Y, Ds, # angle,
                      optimize_weights) {
    toreplace <- sample(1:length(vs), 1)
    replacewith <- sample(as.numeric(V(g)[-which(V(g) %in% vs)]), 1)
    vs[toreplace] <- replacewith
    return(vs)
}

                                        # I would like to eliminate this function
#' Select a set of sentinel nodes by optimization
#'
#' Using combinatorial simulated annealing, select a set of sentinel nodes. Optionally, optimize node weights.
#'
#' @param n The number of nodes in the node set
#' @param y The system mean state
#' @param Y A matrix, the full system state
#' @param Ds The range of control parameter values (length(Ds) should equal length(y) and nrow(Y))
#' @param maxit Defaults NULL. If default, standard values will be used, depending on the number of nodes. Passed to optim().
#' @param optimize_weights Logical for whether or not to optimize node weights.
#' @param trace Logical, passed to optim(). Not working right now.
#'
#' @return
#' A list of three elements:
#' - "vs" is the sorted (by degree) vector of selected nodes,
#' - "error" is the approximation error of the node set,
#' - and "ks" is the degree sequence.
#' @export
experiment <- function(n, y, Y, Ds, maxit = NULL, # angle = FALSE,
                       optimize_weights = FALSE, trace = FALSE) {

    if(is.null(maxit)) {
        maxit <- switch(
            optimize_weights+1, # convert from {0, 1} to {1, 2}
            50*ncol(Y), # more iterations if not optimizing weights
            25*ncol(Y)) # fewer iterations if optimizing weights
    }

    vs <- sample(as.numeric(V(g)), n)
    result <- optim(
        vs, obj_fn, update_vs, y = y, Y = Y, Ds = Ds,
        # angle = angle,
        optimize_weights = optimize_weights,
        method = "SANN", control = list(trace = trace, maxit = maxit, temp = 10)
    )

    vs <- result$par
    error <- result$value
    ks <- k[vs]

    sorting <- order(ks)
    
    if(optimize_weights) {
        ws <- quadoptm(vs, y, Y)
        return(list(vs = vs[sorting], error = error, ks = ks[sorting], ws = ws[sorting]))
    } else {
        return(list(vs = vs[sorting], error = error, ks = ks[sorting]))
    }
}

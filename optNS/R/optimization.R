weighted.mean <- stats::weighted.mean

## #' @import stats
## library(stats)
#' @import igraph
library(igraph)
#' @import ROI
library(ROI)
#' @import ROI.plugin.qpoases
library(ROI.plugin.qpoases)

### Objective function
## z is the reduction, y is the mean system state
### ^-- both of these are over all values of the bifurcation parameter (i.e., index by l)
## Y is the full matrix of x_{i,l}^*
## This is unfortunately confusing notation, but is used consistently throughout this analysis.
#' Objective function
#'
#' Compute the objective function for optimization: the mean square error between the approximation z and the observed mean y over the L values of the control parameter. 
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

#' Calculate the approximated mean state for a sample of nodes
#'
#' Given a sample of nodes, compute the objective function: the weighted or unweighted mean state of the sample of nodes, normalized appropriately. Can be used on its own or inside optim(). 
#' @param vs A sample of nodes
#' @param y The ground truth mean state
#' @param Y A matrix
#' @param optimize_weights Logical for whether or not to optimize node weights.
#' @param ws A vector of node weights, if already obtained
#' @param ... Container for arguments passed from optim()
#'
#' @details
#' The primary purpose of this function is to handle various features of the provided node set and return the approximation error for that node set. By default this function will 
#'
#' For the vector `vs`, it is safer to use the numeric indices of the graph to avoid confusion with vertex names.
#' 
#' @return A numeric vector of approximated system states at each value of the bifurcation parameter
#' 
#' @seealso [calc_obj()]
#' @export
obj_fn <- function(vs, y, Y, optimize_weights = FALSE, ws = NULL, ...) {
    Z <- as.matrix(Y[, vs])

    if(optimize_weights) {
        if(is.null(ws)) ws <- quadoptm(vs, y, Y)
        z <- apply(Z, 1, weighted.mean, ws)
    } else {
        z <- rowMeans(Z)
    }

    calc_obj(z, y)
}



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
    Q0 <- (t(Yp)%*%Yp) # 2* # I had to include the 2, maybe to match the 2 in a0
    a0 <- -(t(y)%*%Yp) # 2*

    qp <- OP(
        Q_objective(Q0, a0),
        L_constraint(rep(1, n), "==", 1)
    )
    qpsol <- ROI_solve(qp, "qpoases")

    return(solution(qpsol))
}

#' Select a new node
#'
#' Given a starting sample of nodes, select a random node from the node set and replace it with a random node not in the node set. Primarily for use inside optim(). 
#' @param vs The selected nodes (i.e., the columns of Y)
#' @param g The network.
#' @param y The mean state at each level of the bifurcation parameter
#' @param Y The state of each node (in columns) at each level of the bifurcation parameter (rows)
#' @param optimize_weights Logical for whether or not to optimize node weights.
#'
#' @return A numeric vector of node indices
#' @export
update_vs <- function(vs, g, y, Y, optimize_weights) {
    toreplace <- sample.local(1:length(vs), 1)
    replacewith <- sample.local(as.numeric(V(g)[-which(V(g) %in% vs)]), 1)
    vs[toreplace] <- replacewith
    return(vs)
}

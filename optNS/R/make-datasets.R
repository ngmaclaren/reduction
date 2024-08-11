#' Make a list of node sets
#'
#' Make a list of node sets (and some additional information) according to a specified algorithm.
#'
#' @param ntrials The number of independent realizations of the node set type
#' @param ns.type One of 'opt', 'rand', 'fixed', 'constr', 'quant', or 'comm'
#' @param ncores Passed to mclapply()
#' @param ... Arguments to be passed to the selector function.
#'
#' @details
#' This function generates several realizations of a given algorithm, optionally with optimized node weights. See the manual page for each function for the required arguments.
#'
#' Algorithms select nodes by the following methods: 'opt', combinatorial simulated annealing; 'rand', uniformly at random; 'fixed', uniformly at random, except that a provided vector of nodes specifies the degree sequence; 'constr', uniformly random, except that the largest 5\% of nodes, by degree, are removed from consideration; 'quant', as constr, but additionally selecting one node at random from n 'bins' of the degree distribution; 'comm', as constr, but preferring to select nodes from different communities.
#' @return A list of node sets.
#' @seealso The family of functions which select node sets, which can be found here: [select_optimized()].
#' @export
make_dataset <- function(ntrials,
                         ns.type = c("opt", "rand", "fixed", "constr", "quant", "comm", "knnconstr"),
                         ncores = 1,
                         ...) {
    ns.type <- match.arg(ns.type)

    selector <- switch(
        ns.type,
        opt = select_optimized,
        rand = select_randomized,
        fixed = select_bydegseq,
        constr = select_constrained,
        quant = select_quantiled,
        comm = select_bycommunity,
        knnconstr = select_knn_constrained
    )

    with(list(args = list(...)), {
        mclapply(seq_len(ntrials), function(x) do.call(selector, args), mc.cores = ncores)
    })
}

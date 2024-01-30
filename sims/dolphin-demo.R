                                        # Options
library(optparse)
optionlist <- list(
    make_option(
        c("-v", "--verbose"), action = "store_true", default = FALSE,
        help = "Show optimization details. Setting this flag sets `trace = TRUE` in `optim()` and `mc.silent = FALSE` in `mclapply()`."
    ),
    make_option(
        c("-s", "--save-plots"), action = "store_true", default = FALSE,
        help = "Store output figures as PDFs. If tabular output is created, saves to a text file. If FALSE [default], will display plot in an R graphics device (and silently fail to do so if graphics devices are not available.) Ignored for simulation files."
    ),
    make_option(
        c("-w", "--optimize-weights"), action = "store_true", default = FALSE,
        help = "Use quadratic programming to optimize weights of nodes selected by simulated annealing. Substantially increases run time (approximately 4--5x)."
    )
)

args <- parse_args(
    OptionParser(option_list = optionlist),
    convert_hyphens_to_underscores = TRUE
)

verbose <- args$verbose
save_plots <- args$save_plots
optimize_weights <- args$optimize_weights
outputfile <- switch(
    optimize_weights + 1,
    "../data/dolphin-demo.RData", # FALSE
    "../data/dolphin-demo-weighted.RData" # TRUE
)

                                        # Parallel processing, including random number set up
library(parallel)
ncores <- detectCores() - 1
RNGkind("L'Ecuyer-CMRG")
set.seed(12345)
                                        # analysis libraries
library(igraph)
library(deSolve)
                                        # Local analysis functions
source("../src/functions.R")


                                        # Load data, set up variables
load("../data/dolphin.rda")
load("../data/fullstate-dolphin.rda")

g <- dolphin
N <- vcount(g)
A <- as_adj(g, "both", sparse = FALSE)
k <- degree(g)

Y <- fullstate$dw
y <- rowMeans(Y)
ns <- 1:floor(log(N))

                                        # GBB and DART
GBB.b <- mean(k^2)/mean(k)
GBB <- with(doublewell_parms, {
    solve_doublewell(x = x.init, r = r, Ds = Ds, A = matrix(1), b = GBB.b)
})
GBB.obs <- apply(Y, 1, function(x) mean(k*x)/mean(k))

eigs <- eigen(A, symmetric = TRUE)
DART.a <- eigs$vectors[, 1]/sum(eigs$vectors[, 1])
DART.b <- eigs$values[1]
DART <- with(doublewell_parms, {
    solve_doublewell(x = x.init, r = r, Ds = Ds, A = matrix(1), b = DART.b)
})
DART.obs <- apply(Y, 1, function(x) sum(DART.a*x))

                                        # Node set selection
solns <- lapply(ns, function(n) {
    with(doublewell_parms, {
        experiment(n, y, Y, Ds, optimize_weights = optimize_weights, trace = verbose)
    })
})

                                        # For n=1 and n=2, the exact solution
find_exact <- function(n, ...) {
    VS <- combn(N, n)
    errors <- apply(VS, 2, obj_fn, y, Y, Ds, optimize_weights = optimize_weights)

    idx <- which.min(errors)
    vs <- VS[, idx]
    error <- errors[idx]
    ks <- k[vs]

    if(optimize_weights) ws <- quadoptm(vs, y, Y)

    if(optimize_weights) {
        return(list(vs = vs, error = error, ks = ks, ws = ws))
    } else {
        return(list(vs = vs, error = error, ks = ks))
    }
}
exact <- list(find_exact(1, y, Y, Ds), find_exact(2, y, Y, Ds))

save.image(outputfile)

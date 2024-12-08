This repository contains code to select optimized sentinel node sets which can closely approximate the state of all network nodes, as well as simulation and analysis code which produced the findings in our manuscript. See below for an example analysis.

## optNS 

The `optNS` package is a small R package which implements our optimization algorithm; the name of the package is simply short for "optimize node sets." The `./optNS/` directory contains the code, which you can build with
```sh
R CMD build optNS
```
The built package is also included as `optNS_0.1-1.tar.gz`. To install, use
```sh
R CMD INSTALL optNS_0.1-1.tar.gz
```

Please contact me with any questions.

## Shell scripts for high-performance computing clusters

Our method should run in reasonable time for small networks on a standard laptop (I use a laptop with four Intel i3-5010U CPUs at 2.10GHz). In order to support larger networks we used a high-performance computing cluster to simulate node states across ranges of control parameters (our "ground truth" simulations) and to optimize node sets. The shell scripts we used to do this are in `./shell/` and will need to be modified for your parallel computing environment.

## Dependencies

In addition to the `optNS` package, we use an in-house differential equation simulation package called [sdn](https://github.com/ngmaclaren/sdn). We use the following standard R repository packages:

- [deSolve](https://cran.r-project.org/package=deSolve)
- [ROI](https://cran.r-project.org/package=ROI)
- [ROI.plugin.qpoases](https://cran.r-project.org/package=ROI.plugin.qpoases)
- [igraph](https://cran.r-project.org/package=igraph)
- [parallel](https://cran.r-project.org/doc/manuals/r-release/fullrefman.pdf)
- [optparse](https://cran.r-project.org/package=optparse)
- [latex2exp](https://cran.r-project.org/package=latex2exp) (for the plots in `./analysis`)

## Example

What follows is a simple demonstration of our method using the functions in `optNS`. The code below is also included in the repository as `demo.R`.

As tested, we rely on a minimum of the following packages (setting the number of cores, color palette, and random seed for convenience):
```R
library(parallel)
ncores <- detectCores()-1
library(igraph)
library(optNS)
palette("R4")
set.seed(1234)
```
We have already simulated $x_i^*$ values for the various networks we use in the manuscript, so for convenience let's use one of those. We'll use the coupled double-well dynamics on the dolphin social network. To run our algorithm, we need the number of nodes in the node set `n`, the network `g`, and the "ground truth" simulation results to use for optimization. The matrix `Y` below is the "ground truth" simulation results for the coupled double-well dynamics on the dolphin network: the rows correspond to control parameter values (here, $D \in [0, 1]$) and the columns to the nodes in the dolphin network. 
```R
network <- "dolphin"
dynamics <- "doublewell"

fullstatefile <- paste0("./data/fullstate-", network, ".rds") # the stored "ground truth" simulation

g <- readRDS(paste0("./data/", network, ".rds")) # the dolphin network
N <- vcount(g)

Y <- readRDS(fullstatefile)[[dynamics]] ## the coupled double-well "ground truth"
y <- rowMeans(Y)
n <- floor(log(N))
```
To produce one optimized nodeset we call `optNS::select_optimized()`
```R
system.time(
    soln <- select_optimized(n, g, y, Y)
)
```
which took about 2.4 seconds on my laptop. Note that larger networks require far more time and/or computing resources: we run the simulated annealing algorithm for $50N$ iterations by default. 

The output of `select_optimized()` is a list of vertex indices, approximation error, and degree sequence; optimized node weights are null unless you ask `select_optimized()` to use the quadratic optimization procedure.
```R
str(soln)
```

To compute the approximation we select the columns from `Y` that correspond to the nodes our optimization algorithm returned:
```R
Z <- Y[, soln$vs]
z <- rowMeans(Z)
```

We can visualize the results as follows:
```R
Ds <- seq(0, 1, length.out = 100)

matplot(
    Ds, Y,
    type = "l", lty = 1, lwd = 0.25, col = "gray50",
    xlab = "D", ylab = expression(x[i]), font.lab = 3
)
lines(Ds, y, lty = 1, lwd = 8, col = 1)
matlines(Ds, Z, lty = 1, lwd = 4, col = 2)
lines(Ds, z, lty = 1, lwd = 8, col = 3)
```

To make several node sets, we need a bit more time and a standardized function `optNS::make_dataset()` which can produce optimized or random node sets, among other options.
```R
ntrials <- 10

system.time(
    opts <- make_dataset(
        ntrials = ntrials, ns.type = "opt", ncores = ncores,
        n = n, g = g, y = y, Y = Y
    )
)
```
That took about 12.5 seconds on my laptop. Making random node sets takes very little time.
```R
rands <- make_dataset(
    ntrials = 10*ntrials, ns.type = "rand", ncores = ncores,
    n = n, g = g, y = y, Y = Y
)
```
`optNS` also includes some convenience functions to retrieve information from lists of node sets. We can use those functions to make a plot that summarizes our simulation results:
```R
error <- list(
    opt = get_error(opts),
    rand = get_error(rands)
)

xlim <- range(unlist(error))

par(mar = c(4, 8, 1, 1))
plot(NULL, xlim = xlim, ylim = c(0.5, 2.5), xlab = "Approximation error", ylab = "", log = "x", yaxt = "n")
points(error$opt, jitter(rep(2, length(error$opt)), amount = 0.1), pch = 1, col = 2)
points(error$rand, jitter(rep(1, length(error$rand)), amount = 0.1), pch = 0, col = 1)
axis(2, at = c(2, 1), tick = FALSE, labels = c("Optimized", "Random"), las = 2)
```

library(parallel)
ncores <- detectCores()-1
library(igraph)
library(optNS)

network <- "dolphin" # "email"
dynamics <- "doublewell"

fullstatefile <- paste0("../data/fullstate-", network, ".rds")

g <- readRDS(paste0("../data/", network, ".rds"))
N <- vcount(g)
AL <- as_adj_list(g, "all")

Y <- readRDS(fullstatefile)[[dynamics]]
y <- rowMeans(Y)

tc_ns <- 1:30
ntrials <- 5

tc_runtimes <- matrix(NA, nrow = length(tc_ns), ncol = ntrials)
tc_results <- list()

tc_sim <- function(n) {
    runtime <- system.time(
        result <- select_optimized(n, g, y, Y)
    )
    c(runtime = runtime[3], n = n)
}

do.call(rbind, mclapply(seq(3), function(i) tc_sim(3)))

for(i in seq_along(tc_ns)) {
    for(j in seq(ntrials)) {
        n <- tc_ns[i]
        runtime <- system.time(
            result <- select_optimized(n, g, y, Y)
        )
        tc_runtimes[i, j] <- runtime[3]
        tc_results[[i]] <- result
    }
}
rm(i, j, n, runtime)


## Now, vary maxit
maxits <- seq(10, 100, by = 10)*N
hm_ns <- seq(5, 30, by = 5)
hm_results <- list()

for(i in seq_along(maxits)) {
    maxit <- maxits[i]
    hm_results[[i]] <- list()
    for(j in seq_along(hm_ns)) {
        n <- hm_ns[j]
        hm_results[[i]][[j]] <- make_dataset(
            ntrials = ntrials, ns.type = "opt", ncores = ncores,
            n = n, g = g, y = y, Y = Y, maxit = maxit
        )
    }
}
rm(i, j, n, maxit)

save.image(paste0("tcv-", network, ".RData"))

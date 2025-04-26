library(parallel)
ncores <- detectCores()-2
library(igraph)
library(sdn)
library(optNS)

datadir <- "/projects/academic/naokimas/neil/brains-ns300/"
datafiles <- list.files(path = datadir, pattern = ".txt")

                                        # for testing
## datafile <- datafiles[10]
## datafile <- "./test/101006-50.txt"

make_nodesets <- function(datafile, use.A2 = FALSE) {
    calc_obj <- function(z, y) sum((z - y)^2)/length(y)
    obj_fn <- function (vs, y, Y, optimize_weights = FALSE, ws = NULL, ...) {
        Z <- as.matrix(Y[, vs])
        if (optimize_weights) {
            if (is.null(ws)) 
                ws <- quadoptm(vs, y, Y)
            z <- apply(Z, 1, weighted.mean, ws)
        }
        else {
            z <- rowMeans(Z)
        }
        calc_obj(z, y)
    }

    df <- read.table(paste0(datadir, datafile), sep = " ")
    nr <- nrow(df)
    pos <- which(datafiles %in% datafile)

    X.train <- as.matrix(df[seq_len(nr-100), ])
    X.test <- as.matrix(df[seq(nr-100+1, nr, by = 1), ])

    avgxi <- colMeans(X.test)
    mu <- mean(avgxi) # mean(X.test)
    sigma <- sd(avgxi) # sd(X.test)
    threshold <- mu + 5*sigma # 4
    test <- any(abs(avgxi) >= threshold) # colMeans(X.train)
    if(isFALSE(test)) flag <- "use" else flag <- "dontuse"
    ##return(flag)
##}
    if(use.A2) {
        A <- readRDS("A2.rds")[[pos]]
    } else {
        Cmat <- cor(X.train, method = "pearson")
        A <- Cmat
    }
    A[which(A < 0, arr.ind = TRUE)] <- 0
    if(use.A2) A <- A/max(A)

    g <- graph_from_adjacency_matrix(A, "undirected", weighted = TRUE, diag = FALSE)
    return(g)
}


##     N <- vcount(g)
##     n <- floor(log(N))

##     times <- 0:15
##     control <- list(times = times, ncores = ncores)
##     params <- c(.doublewell, list(A = A))
##     params$use.matrix <- TRUE
##     Ds <- seq(0, 1, length.out = 100)
##     Y <- solve_in_range(Ds, "D", doublewell, rep(params$xinit.low, N), params, control, "ode")
##     y <- rowMeans(Y)

##     opts <- make_dataset(
##         ntrials = 100, ns.type = "opt", ncores = ncores, n = n, g = g, y = y, Y = Y
##     )
##     rands <- make_dataset(
##         ntrials = 100, ns.type = "rand", ncores = ncores, n = n, g = g, y = y, Y = Y
##     )
##     testerror <- list(
##         opt = sapply(opts, function(ns) obj_fn(ns$vs, rowMeans(X.test), X.test)),
##         rand = sapply(rands, function(ns) obj_fn(ns$vs, rowMeans(X.test), X.test))
##     )
##     t.test(testerror$opt, testerror$rand, "less")
    
##     return(list(network = g, opts = opts, rands = rands, testerror = testerror, flag = flag))
## }


## datalist <- lapply(datafiles, make_nodesets)
## saveRDS(datalist, "datalist-large.rds")

## check <- unlist(mclapply(datafiles, make_nodesets, mc.cores = ncores))
## saveRDS(check, "check5.rds")

graphlist <- lapply(datafiles, make_nodesets)
saveRDS(graphlist, "graphlist-large.rds")

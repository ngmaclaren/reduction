library(parallel)
ncores <- detectCores()-1
library(igraph)
library(optNS)
KL <- philentropy::KL

dolphin.env <- new.env()
ba.env <- new.env()

with(dolphin.env, {
    g <- readRDS("../data/dolphin.rds")
    N <- vcount(g)
    k <- degree(g)
    Y <- readRDS("../data/fullstate-dolphin.rds")$doublewell
    y <- rowMeans(Y)
    ns <- 1:12
    ntrials <- 100

    opts <- readRDS("../data/dolphin-1-to-12.rds")
    allrands <- replicate(
        ntrials,
        lapply(ns, function(n) {
            make_dataset(
                ntrials = ntrials, ns.type = "rand", ncores = ncores,
                n = n, g = g, y = y, Y = Y
            )
        }),
        FALSE
    )

    kposs <- seq_len(max(k))
    counts.orig <- sapply(kposs, function(ki) sum(k == ki))
    counts.opts <- lapply(opts, function(dl) sapply(kposs, function(ki) sum(as.numeric(get_ks(dl)) == ki)))
    KLs.opts <- sapply(counts.opts, function(opt) KL(rbind(opt, counts.orig), est.prob = "empirical"))
    errors.opts <- sapply(opts, get_error) # n in cols, ntrials in rows

    counts.rands <- lapply(
        allrands,
        function(rands) lapply(rands, function(dl) sapply(kposs, function(ki) sum(as.numeric(get_ks(dl)) == ki)))
    )

    KLs.rands <- lapply(
        counts.rands,
        function(counts) sapply(counts, function(rand) KL(rbind(rand, counts.orig), est.prob = "empirical"))
    )
    KLs.rands <- do.call(rbind, KLs.rands)

    errors.rands <- lapply(allrands, function(rands) sapply(rands, function(dl) mean(get_error(dl))))
    errors.rands <- do.call(rbind, errors.rands)
})

with(ba.env, {
    g <- readRDS("../data/ba.rds")
    N <- vcount(g)
    k <- degree(g)
    Y <- readRDS("../data/fullstate-ba.rds")$doublewell
    y <- rowMeans(Y)
    ns <- 1:12
    ntrials <- 100

    opts <- readRDS("../data/ba-1-to-12.rds")
    allrands <- replicate(
        ntrials,
        lapply(ns, function(n) {
            make_dataset(
                ntrials = ntrials, ns.type = "rand", ncores = ncores,
                n = n, g = g, y = y, Y = Y
            )
        }),
        FALSE
    )

    kposs <- seq_len(max(k))
    counts.orig <- sapply(kposs, function(ki) sum(k == ki))
    counts.opts <- lapply(opts, function(dl) sapply(kposs, function(ki) sum(as.numeric(get_ks(dl)) == ki)))
    KLs.opts <- sapply(counts.opts, function(opt) KL(rbind(opt, counts.orig), est.prob = "empirical"))
    errors.opts <- sapply(opts, get_error) # n in cols, ntrials in rows

    counts.rands <- lapply(
        allrands,
        function(rands) lapply(rands, function(dl) sapply(kposs, function(ki) sum(as.numeric(get_ks(dl)) == ki)))
    )

    KLs.rands <- lapply(
        counts.rands,
        function(counts) sapply(counts, function(rand) KL(rbind(rand, counts.orig), est.prob = "empirical"))
    )
    KLs.rands <- do.call(rbind, KLs.rands)

    errors.rands <- lapply(allrands, function(rands) sapply(rands, function(dl) mean(get_error(dl))))
    errors.rands <- do.call(rbind, errors.rands)
})

save.image("kld-fig.RData")

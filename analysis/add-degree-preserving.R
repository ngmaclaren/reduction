library(parallel)
ncores <- detectCores()-2
library(igraph)
library(sdn)
library(optNS)

flag <- readRDS("check5.rds")
graphlist <- readRDS("./graphlist-large.rds")[which(flag == "use")]
## fullstates <- readRDS("./fullstates-fmri.rds")[which(flag == "use")]
dl <- readRDS("datalist-large.rds")[which(flag == "use")]
datadir <- "/projects/academic/naokimas/neil/brains-ns300/"
datafiles <- list.files(path = datadir, pattern = ".txt")[which(flag == "use")]

make_fixeds <- function(datafile, opts, g) {
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

    N <- vcount(g)
    n <- floor(log(N))

    df <- read.table(paste0(datadir, datafile), sep = " ")
    nr <- nrow(df)
    pos <- which(datafiles %in% datafile)

    X.train <- as.matrix(df[seq_len(nr-100), ])
    X.test <- as.matrix(df[seq(nr-100+1, nr, by = 1), ])

    ## Cmat <- cor(X.train, method = "pearson")
    ## A <- Cmat
    ## A[which(A < 0, arr.ind = TRUE)] <- 0
    ## diag(A) <- 0

    best <- opts[[which.min(get_error(opts))]]

    fixeds <- make_dataset(
        ntrials = 100, ns.type = "fixed", ncores = ncores,
        n = n, g = g, comps = best$vs, y = rowMeans(X.test), Y = X.test
    )

    testerror <- sapply(fixeds, function(ns) obj_fn(ns$vs, rowMeans(X.test), X.test))

    return(list(fixeds = fixeds, testerror = testerror))
}

##test <- make_fixeds(datafiles[1], dl[[1]]$opts, graphlist[[1]])

## datafile, opts, g
opts <- lapply(dl, `[[`, "opts")
newdata <- mapply(make_fixeds, datafiles, opts, graphlist)
    
## dl$fixeds <- lapply(newdata, `[[`, "fixeds")
## dl$testerror$fixed <- newdata$testerror

saveRDS(newdata, "datalist-large-dp.rds")
    
## then add fixed to testerror as well

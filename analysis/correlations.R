## need to adjust this file manually for the different networks. Options are celegans and proximity for now.

if(getwd() != "/user/neilmacl/Documents/reduction/analysis") if(interactive()) setwd("./analysis")

mut.env <- new.env()
dw.env <- new.env()
sis.env <- new.env()

                                        # Load the network into the top env
load("../data/celegans.rda")
                                        # load the simulation spaces into their own envs
load("../data/sentinel-nodesets-celegans-mutualistic.RData", mut.env)
load("../data/sentinel-nodesets-celegans-dw.RData", dw.env)
load("../data/sentinel-nodesets-celegans-SIS.RData", sis.env)

output_data <- function(env) {
    with(env, {
        optsnodes <- as.numeric(sapply(opts, `[[`, "vs"))
        ## count the number of times that node i is contained in an optimized node set
        counts <- sapply(1:N, function(i) sum(optsnodes == i))
        p <- counts/ntrials
        ## calculate the mean error for random node sets that contain a node
        e <- sapply(1:N, function(i) {
            whichsets <- sapply(rands, function(rand) i %in% rand$vs) # TRUE or FALSE
            if(sum(whichsets) > 0) {
                mean(sapply(rands[whichsets], `[[`, "error"))
            } else {
                NA
            }
        })  
    })

    return(list(p = with(env, p), e = with(env, e)))
}

mut.data <- output_data(mut.env)
dw.data <- output_data(dw.env)
sis.data <- output_data(sis.env)

save(mut.data, dw.data, sis.data, celegans, file = "celegans-data.rda")

        
## with(A.env, optsnodes <- as.numeric(sapply(opts, `[[`, "vs")))
## with(B.env, optsnodes <- as.numeric(sapply(opts, `[[`, "vs")))

## with(A.env, {
##     counts <- sapply(1:N, function(i) sum(optsnodes == i))
##     p <- counts/ntrials
## })
## with(B.env, {
##     counts <- sapply(1:N, function(i) sum(optsnodes == i))
##     p <- counts/ntrials
## })

## cor(with(A.env, p), with(B.env, p))

## pdf("plotp.pdf")
## plot(with(A.env, counts), with(B.env, counts),
##      xlab = "Mutualistic dynamics", ylab = "SIS dynamics")
## dev.off()

## ## what is the average approximation error for node i when i is included in a random node set?
## with(A.env, {
##     e <- sapply(1:N, function(i) {
##         whichsets <- sapply(rands, function(rand) i %in% rand$vs) # TRUE or FALSE
##         if(sum(whichsets) > 0) {
##             mean(sapply(rands[whichsets], `[[`, "error"))
##         } else {
##             NA
##         }
##     })
## })
## with(B.env, {
##     e <- sapply(1:N, function(i) {
##         whichsets <- sapply(rands, function(rand) i %in% rand$vs) # TRUE or FALSE
##         if(sum(whichsets) > 0) {
##             mean(sapply(rands[whichsets], `[[`, "error"))
##         } else {
##             NA
##         }
##     })
## })

## cor(with(A.env, e), with(B.env, e), use = "pairwise.complete.obs")

## pdf("plote.pdf")
## plot(with(A.env, e), with(B.env, e), # log = "xy",
##      xlab = "Mutualistic dynamics", ylab = "SIS dynamics")
## dev.off()




## ## Do some regression analysis here.
## ## degree, local clustering, knn
## library(igraph)
## load("../data/proximity.rda")
## g <- proximity

## k <- degree(g)
## bc <- betweenness(g, directed = FALSE)
## knn <- knn(g)$knn
## lcl <- transitivity(g, "localundirected")

## Amod <- glm(with(A.env, e) ~ k + knn + lcl, family = Gamma("log"))
## Bmod <- glm(with(B.env, e) ~ k + knn + lcl, family = Gamma("log"))



## A.e <- with(A.env, e)
## B.e <- with(B.env, e)
## save(A.e, B.e, k, bc, knn, lcl, g, file = "proximity-data.rda")

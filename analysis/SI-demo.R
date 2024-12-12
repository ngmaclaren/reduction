                                        # packages #
library(parallel)
ncores <- detectCores()-1
library(igraph)
library(sdn)
library(optNS)

                                        # functions #
                                        # Chebyshev/GBB, including coefs
polyval <- pracma::polyval
tfun <- function(x, xmax) (2*x/xmax) - 1
chebyG <- function(x, order, dynamics, network) {
    xmax <- cheby[[dynamics]][[network]]$xmax
    coef <- cheby[[dynamics]][[network]][[paste0("order", order)]]
    polyval(coef, tfun(x, xmax))
}
mutualisticGBB <- function(t, x, params) {
    with(params, {
        dx <- ifelse(
            x > 0,
            B + x*(1 - (x/K)) * ((x/C) - 1) + D*b*(x^2/(Dtilde + (E + H)*x)),
            0
        )
        return(list(c(dx)))
    })
}
mutualisticCheby <- function(t, x, params) {# order, dynamics, network
    with(params, {
        dx <- ifelse(
            x > 0,
            B + x * (1 - (x/K)) * ((x/C) - 1) + D*b*chebyG(x, order, dynamics, network),
            0
        )
        return(list(c(dx)))
    })
}
generegCheby <- function(t, x, params) {
    with(params, {
        dx <- ifelse(
            x > 0,
            -B*(x^f) + D*b*chebyG(x, order, dynamics, network),
            0
        )
        return(list(c(dx)))
    })
}
                                        # Cheby coefs.
cheby <- list( 
    mutualistic = list(
        hk = list(
            xmax = 26.870522,
            order3 = c(t6 = 0.018441, t5= -0.539394, t4 = 2.570592, t3 = -0.184549, t2 = -1.193048, t1 = 12.199392,
                       t0 = 10.025042),
            order2 = c(t4 = 0.387120, t3 = -3.491421, t2 = 0.486049, t1 = 14.376137, t0 = 10.025042)
        ),
        dolphin = list(
            xmax = 13.80085,
            order3 = c(t6 = 0.003160, t5 = -0.120503, t4 = 0.753642, t3 = -0.432709, t2 = 0.281885,
                       t1 = 5.652725, t0 = 4.051398),
            order2 = c(t4 = 0.120214, t3 = -1.305299, t2 = 0.764067, t1 = 6.239385, t0 = 4.051398)
        ),
        ba = list(
            xmax = 25.14951,
            order3 = c(t6 = 0.015810, t5 = -0.472493, t4 = 2.302986, t3 = -0.258229, t2 = -0.947378,
                       t1 = 11.349985, t0 = 9.201113),
            order2 = c(t4 = 0.349292, t3 = -3.191047, t2 = 0.553464, t1 = 13.283822, t0 = 9.201113)
        )
    ),
    genereg = list(
        drosophila = list(
            xmax = 124.5753,
            order3 = c(t3 = 0.450857, t2 = -0.240146, t1 = -0.212957, t0 = 1.056467),
            order2 = c(t2 = -0.240146, t1 = 0.125186, t0 = 1.056467)
        ),
        dolphin = list(
            xmax = 11.60787,
            order3 = c(t3 = 0.586945, t2 = -0.503663, t1 = -0.076728, t0 = 1.035964),
            order2 = c(t2 = -0.503663, t1 = 0.363481, t0 = 1.035964)
        ),
        ba = list(
            xmax = 72.72385,
            order3 = c(t3 = 0.547207, t2 = -0.303240, t1 = -0.247971, t0 = 1.068138),
            order2 = c(t2 = -0.303240, t1 = 0.162434, t0 = 1.068138)
        )
    )
)

                                        # DART
## Nomenclature:
## my x = DART's R (1D state variable)
## my b = DART's alpha (ith eigenvalue)
## my beta = DART's beta (structural parameter, multiplied by x to replace x_i)
## my a = DART's a (ith eigenvector)
calc_beta <- function(i, g, alpha = NULL, a = NULL) {
    if(is.null(alpha) & is.null(a)) {
        A <- as_adjacency_matrix(g, "both", sparse = FALSE)
        eigs <- eigen(A, symmetric = isSymmetric(A))
        alpha <- eigs$values[i]
        a <- eigs$vectors[, i]/sum(eigs$vectors[, i])
    }
    K <- diag(degree(g))
    return(as.numeric((1/alpha)*((t(a) %*% K %*% a)/(t(a) %*% a))))
}
mutualistic1D <- function(t, x, params) {
    with(params, {
        coupling <- D*b*((beta*x*x))/(Dtilde + E*beta*x + H*x)
        dx <- ifelse(x > 0, B + x * (1 - (x/K)) * ((x/C) - 1) + coupling + u, 0)
        return(list(c(dx)))
    })
}
genereg1D <- function(t, x, params) {
    with(params, {
        dx <- ifelse(x > 0, -B*(x^f) + D*b*(x^h/(1 + x^h)) + u, 0)
        return(list(c(dx)))
    })
}

                                        # plotting
plot_ns <- function(Ds, Y, ...) {
    matplot(
        Ds, Y, type = "l", lty = 1, lwd = 0.5, col = "#babdb6",
        xlab = "D", ylab = "x", cex.lab = labelsize, cex.axis = labelsize, ...
    )
}
labelsize <- 1.75

                                        # selectors (dynamics, network, ns size) #
                                        # Can make this accept optparse here.
dynamics <- "genereg" # "mutualistic"
network <- "drosophila" # "hk"
params <- get(paste0(".", dynamics))
xinit <- switch(dynamics, mutualistic = params$xinit.low, genereg = params$xinit.high)

outfile <- paste0("SI-comps-", dynamics, ".RData")

                                        # set up the network, other needed variables #
g <- readRDS(paste0("../data/", network, ".rds"))
N <- vcount(g)
A <- as_adjacency_matrix(g, "both", sparse = FALSE)
k <- degree(g)
Y <- readRDS(paste0("../data/fullstate-", network, ".rds"))[[dynamics]]
y <- rowMeans(Y)
eigs <- eigen(A, TRUE)
a <- sapply(1:3, function(i) eigs$vectors[, i]/sum(eigs$vectors[, i])) # check if needs sign switch
DART.b <- sapply(1:3, function(i) eigs$values[i])
beta <- sapply(1:3, function(i) calc_beta(i, g, alpha = DART.b[i], a = a[, i]))
GBB.b <- mean(k^2)/mean(k)
control <- list(times = 0:15, ncores = ncores)

### Our method ###

ns <- 1:4
system.time(solns <- mclapply(ns, function(n) select_optimized(n, g, y, Y), mc.cores = ncores))
get_error(solns)


### GBB/Chebyshev ###

                                        # Check Cheby
x <- seq(0, cheby[[dynamics]][[network]]$xmax, length.out = 1000)
fx <- switch(
    dynamics,
    genereg = x^2/(1 + x^2),
    mutualistic = with(.mutualistic, x^2/(Dtilde + x))
)
o2 <- chebyG(x, 2, dynamics, network)
o3 <- chebyG(x, 3, dynamics, network)

plot(x, fx, type = "l", col = "black", ylim = range(c(fx, o2, o3)))
lines(x, o2, col = "red")
lines(x, o3, col = "blue")
legend("bottomright", bty = "n", legend = c(2, 3), title = "Order", lty = 1, col = c("red", "blue"))

## ylim <- switch(dynamics, mutualistic = c(-0.2, 0.2))
## xlim <- switch(dynamics, mutualistic = c(0, 1))
## plot(x, fx, type = "l", col = "black", xlim = xlim, ylim = ylim, yaxs = "i", xaxs = "i")
## lines(x, o2, col = "red")
## lines(x, o3, col = "blue")
## marker <- which(o3 > 0)[1]
## abline(v = x[marker], lwd = 0.5, col = "blue")
## ##points(x[marker], o3[marker], pch = 1, col = "blue", cex = 3, lwd = 2)
## abline(h = 0, col = "gray50", lty = 3)
## legend("bottomright", bty = "n", legend = c(2, 3), title = "Order", lty = 1, col = c("red", "blue"))


                                        # theory/reduction
                                        # 1-D
GBB <- solve_in_range(
    params$Ds, "D", switch(dynamics, mutualistic = mutualisticGBB, genereg = genereg1D), xinit,
    params = c(params, list(b = GBB.b)), control = list(times = 0:15, ncores = 1)
)

                                        # 2-D
Cheby2 <- solve_in_range(
    params$Ds, "D", switch(dynamics, mutualistic = mutualisticCheby, genereg = generegCheby), xinit,
    params = c(params, list(b = GBB.b, order = 2, dynamics = "mutualistic", network = "hk")),
    control = control
)

                                        # 3-D
Cheby3 <- solve_in_range(
    params$Ds, "D", switch(dynamics, mutualistic = mutualisticCheby, genereg = generegCheby),
    switch(dynamics, mutualistic = 0.5, genereg = xinit), # see diagnostics above
    params = c(params, list(b = GBB.b, order = 3, dynamics = "mutualistic", network = "hk")),
    control = control, method = "adams", maxsteps = 20000
)

                                        # ground truth/numerical/simulated
GBB.obs <- apply(Y, 1, function(x) mean(k*x)/mean(k))

                                        # error against mean
calc_obj(GBB, GBB.obs)
calc_obj(GBB, y)

calc_obj(Cheby2, GBB.obs)
calc_obj(Cheby2, y)

calc_obj(Cheby3, GBB.obs)
calc_obj(Cheby3, y)
                                        # error against GBB observable

### DART ###

                                        # theory/reduction
DARTs <- sapply(1:3, function(i) {
    a <- a[, i]
    b <- DART.b[i]
    beta <- beta[i]
    solve_in_range(
        params$Ds, "D", switch(dynamics, mutualistic = mutualistic1D, genereg = genereg1D), xinit,
        params = c(params, list(a = a, b = b, beta = beta)), control = control
    )
})

                                        # ground truth/numerical/simulated
DART.obs.parts <- sapply(1:3, function(i) {
    apply(Y, 1, function(x) sum(a[, i]*x))
})
DART.obs <- list(DART.obs.parts[, 1])
DART.obs[[2]] <- rowMeans(DART.obs.parts[, 1:2])
DART.obs[[3]] <- rowMeans(DART.obs.parts[, 1:3])
DART.obs <- do.call(cbind, DART.obs)
                                        # error against mean
apply(DARTs, 2, calc_obj, y = y)

                                        # error against DART observable
sapply(1:3, function(i) calc_obj(DARTs[, i], DART.obs[, i]))

## save image or make plot
save.image(outfile)


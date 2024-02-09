library(ROI)#, lib.loc = "/user/neilmacl/rlocal/")
library(ROI.plugin.qpoases)#, lib.loc = "/user/neilmacl/rlocal/")

lout <- 100 # `L` in the manuscript: the number of samples of the control parameter

### Convert these simulation/solving functions to newer, clearer form? Requires re-running from the beginning, so if need to do that anyway then go ahead.

### Mutualistic Species Dynamics
mutualistic <- function(t, x, parms) {
    with(parms, {
        coupling <- D*rowSums(A*(outer(x, x)/(D.tilde + outer(E*x, H*x, `+`)))) # col
        dx <- B + x*(1 - (x/K))*((x/C) - 1) + coupling
        return(list(c(dx)))
    })
}

solve_mutualistic <- function(x, B, K, C, Ds, D.tilde, E, H, A, times = 0:15) {
    if("parallel" %in% sessionInfo()$basePkgs) {
        ncores <- detectCores() - 1
        results <- mclapply(
            Ds, function(D) {
                ode(x, times, mutualistic,
                    parms = list(B = B, K = K, C = C, D = D, D.tilde = D.tilde, E = E, H = H, A = A))
            }, mc.cores = ncores, mc.silent = TRUE
        )
    } else {
        results <- lapply(
            Ds, function(D) {
                ode(x, times, mutualistic,
                    parms = list(B = B, K = K, C = C, D = D, D.tilde = D.tilde, E = E, H = H, A = A))
            })
    }

    do.call(rbind, lapply(results, function(df) df[nrow(df), -1]))
}

mutualistic_parms <- list(
    x.init = 0.001,
    B = 0.1,
    K = 5,
    C = 1,
    Ds = seq(0, 3, length.out = lout),
    D.tilde = 5,
    E = 0.9,
    H = 0.1,
    times = 0:15
)

### Gene Regulatory Dynamics
genereg <- function(t, x, parms) {
    with(parms, {
        dx <- -B*(x^f) + D*colSums(A*((x^h)/(1 + (x^h))))
        return(list(c(dx)))
    })
}

solve_genereg <- function(x, B, f, h, Ds, A, times = 0:15) {
    if("parallel" %in% sessionInfo()$basePkgs) {
        ncores <- detectCores() - 1
        results <- mclapply(
            Ds, function(D) ode(x, times, genereg, parms = list(B = B, f = f, h = h, D = D, A = A)),
            mc.cores = ncores, mc.silent = TRUE
        )
    } else {
        results <- lapply(Ds, function(D) {
            ode(x, times, genereg, parms = list(B = B, f = f, h = h, D = D, A = A))
        })
    }

    do.call(rbind, lapply(results, function(df) df[nrow(df), -1]))
}

genereg_parms <- list(
    x.init = 2,
    B = 1,
    f = 1,
    h = 2,
    Ds = seq(0, 1, length.out = lout)
)

### Double-well
doublewell <- function(t, x, parms) {
    with(parms, {
        dx <- -(x - r[1])*(x - r[2])*(x - r[3]) + b*D*colSums(A*x)
        return(list(c(dx)))
    })
}

solve_doublewell <- function(x, r, Ds, A, b = 1, times = seq(from = 0, to = 15, by = 1)) {
    if("parallel" %in% sessionInfo()$basePkgs){
        ncores <- detectCores() - 1
        results <- mclapply(
            Ds, function(D) ode(x, times, doublewell, parms = list(r = r, D = D, A = A, b = b)),
            mc.cores = ncores, mc.silent = TRUE
        )
    } else {
        results <- lapply(Ds, function(D) {
            ode(x, times, doublewell, parms = list(r = r, D = D, A = A))
        })
    }
    
    do.call(rbind, lapply(results, function(df) df[nrow(df), -1]))
}

doublewell_parms <- list(
    x.init = 1,
    r = c(1, 3, 5),
    Ds = seq(0, 1, length.out = lout)
)

### SIS
SIS <- function(t, x, parms) {# noiseless version. Needs a check for negative vals if noise
    with(parms, {
        dx <- -μ*x + D*rowSums(A * outer((1 - x), x))
        return(list(c(dx)))
    })
}

solve_SIS <- function(x, μ, Ds, A, times = seq(from = 0, to = 20, by = 1)) {
    if("parallel" %in% sessionInfo()$basePkgs) {
        ncores <- detectCores() - 1
        results <- mclapply(
            Ds, function(D) ode(x, times, SIS, parms = list(μ = μ, D = D, A = A)),
            mc.cores = ncores, mc.silent = TRUE
        )
    } else {
        results <- lapply(Ds, function(D) {
            ode(x, times, SIS, parms = list(μ = μ, D = D, A = A))
        })
    }

    do.call(rbind, lapply(results, function(df) df[nrow(df), -1]))
}

SIS_parms <- list(
    x.init = 0.01,
    μ = 1,
    Ds = seq(0, 1, length.out = lout)
)

### Objective function
## z is the reduction, y is the mean system state
### ^-- both of these are over all values of the bifurcation parameter (i.e., index by l)
## Y is the full matrix of x_{i,l}^*
## This is unfortunately confusing notation, but is used consistently throughout this analysis.

                                        # Monday, December 18, 2023
                                        # Normalize by the grand mean in addition to L
calc_obj <- function(z, y) sum((z - y)^2)/(lout*mean(y))

### Optimization
take_sample <- function(vs, Y) Y[, vs]

                                        # Written this way to support using the Kuramoto model,
                                        # which we decided not to do.
## system_state <- function(x, angle = FALSE) {
##     if(angle) {
##         abs(sum(exp(1i * x))/N)
##     } else {
##         mean(x)
##     }
## }

quadoptm <- function(vs, y, Y) {
    n <- length(vs)
    
    Yp <- Y[, vs]
    Q0 <- 2*(t(Yp)%*%Yp) # I had to include the 2, maybe to match the 2 in a0
    a0 <- -2*(t(y)%*%Yp)

    qp <- OP(
        Q_objective(Q0, a0),
        L_constraint(rep(1, n), "==", 1)
    )
    qpsol <- ROI_solve(qp, "qpoases")

    return(solution(qpsol))
}

reduction.sample <- function(vs, y, Y, Ds, # angle = FALSE,
                             optimize_weights = FALSE) {
    Z <- as.matrix(take_sample(vs, Y))
    if(optimize_weights) {
        w <- quadoptm(vs, y, Y)
        apply(Z, 1, weighted.mean, w)
    } else {
        ##apply(Z, 1, system_state, angle)
        rowMeans(Z)
    }
}

obj_fn <- function(vs, y, Y, Ds, # angle = FALSE,
                   optimize_weights = FALSE) {
    z <- reduction.sample(vs, y, Y, Ds, # angle,
                          optimize_weights)
    calc_obj(z, y)#, Y)
}

update_vs <- function(vs, y, Y, Ds, # angle,
                      optimize_weights) {
    toreplace <- sample(1:length(vs), 1)
    replacewith <- sample(as.numeric(V(g)[-which(V(g) %in% vs)]), 1)
    vs[toreplace] <- replacewith
    return(vs)
}

experiment <- function(n, y, Y, Ds, maxit = NULL, # angle = FALSE,
                       optimize_weights = FALSE, trace = FALSE) {

    if(is.null(maxit)) {
        maxit <- switch(
            optimize_weights+1, # convert from {0, 1} to {1, 2}
            50*ncol(Y), # more iterations if not optimizing weights
            25*ncol(Y)) # fewer iterations if optimizing weights
    }

    vs <- sample(as.numeric(V(g)), n)
    result <- optim(
        vs, obj_fn, update_vs, y = y, Y = Y, Ds = Ds,
        # angle = angle,
        optimize_weights = optimize_weights,
        method = "SANN", control = list(trace = trace, maxit = maxit, temp = 10)
    )

    vs <- result$par
    error <- result$value
    ks <- k[vs]

    sorting <- order(ks)
    
    if(optimize_weights) {
        ws <- quadoptm(vs, y, Y)
        return(list(vs = vs[sorting], error = error, ks = ks[sorting], ws = ws[sorting]))
    } else {
        return(list(vs = vs[sorting], error = error, ks = ks[sorting]))
    }
}

### Monte Carlo functions
select_by_community <- function(n, g, partition, vs = V(g)) {
    stopifnot(is.numeric(n))
    stopifnot(is.igraph(g))
    stopifnot(inherits(partition, "communities"))
    
    stopifnot(length(partition) >= n)

    withheld <- V(g)[-which(V(g) %in% vs)]
    if(length(withheld) > 0) {
        available <- V(g)[-withheld]
    } else {
        available <- vs
    }
    mbr <- membership(partition)
    S <- c()
    m <- c()

    while(length(S) < n) {
        
        v <- sample(available, 1)
        S <- c(S, v)
        m <- c(m, mbr[v])
        
        toremove <- unique(c(which(mbr %in% m), as.numeric(withheld)))
        available <- V(g)[-toremove]
    }

    return(S)
}

select_by_comm_prob <- function(n, g, partition, pvec, vs = V(g)) {
    stopifnot(is.numeric(n))
    stopifnot(is.igraph(g))
    stopifnot(inherits(partition, "communities"))
    stopifnot(length(pvec) == length(partition))

    withheld <- V(g)[-which(V(g) %in% vs)]
    if(length(withheld) > 0) {
        available <- V(g)[-withheld]
    } else {
        available <- vs
    }
    Cs <- seq_len(length(partition))
    mbr <- membership(partition)
    
    C_vec <- sample(Cs, n, TRUE, pvec) # this tells me which communities to take from

    df <- data.frame(v = as.numeric(available), C = as.numeric(mbr[available]))

    sapply(C_vec, function(C) sample(df$v[df$C == C], 1))
    
}


make_dataset <- function(n, # how many nodes are in the node set
                         ntrials, # how many times should we do the experiment
                         bparam, # the bifurcation parameter values (a sequence)
                         y, # make sure this function accepts the right mean state
                         Y, # and full state
                         comps = NULL, # optionally, a degree distribution to fix, provided as vs
                         optimize = FALSE, # should de novo optimized solutions be produced?
                         optimize_weights = FALSE, # should weights be optimized?
                         use_connections = FALSE, # should the ks/knns be used to select nodes?
                         use_quantiles = FALSE, # if use_connections, also split into n by quantile?
                         use_communities = FALSE,
                         verbose = FALSE # for `trace` and `mc.silent`; doesn't seem to do anything
                         ) {
    if(is.numeric(comps)) { # If we have some reference node set to match (e.g. a good node set)
        ks <- k[comps]
        poss <- lapply(ks, function(ki) as.numeric(V(g)[which(k == ki)]))
        VS <- replicate(ntrials, sapply(poss, sample, 1), simplify = FALSE) # fixed-degree
    } else {
        if(optimize) { # fully optimize
            dl <- mclapply(
                1:ntrials, 
                function(x) experiment(
                                n, y, Y, bparam, optimize_weights = optimize_weights,
                                trace = isTRUE(verbose)
                            ),
                mc.cores = ncores, mc.silent = isFALSE(verbose)#FALSE
            )
            return(dl) # can exit with this list of optimized node sets
        } else if(use_connections) { # impose degree/knn based constraints
            library(igraph)
            k <- degree(g)
            knn <- knn(g)$knn

                                        # Effect of knn does not appear to be robust
            ## top10k <- which(k > quantile(k, probs = 0.9))
            ## bottom10knn <- which(knn < quantile(knn, probs = 0.1))
            ## avail <- as.numeric(V(g)[-unique(c(top10k, bottom10knn))])
                                        # Better results with only eliminating the very largest nodes
            top5k <- which(k > quantile(k, probs = 0.95))
            avail <- as.numeric(V(g)[-top5k])
            
            if(use_quantiles) { # impose stratified sampling
                probwidth <- 1/n
                probs <- seq(from = 0, to = 1, by = probwidth)
                quantiles <- quantile(k[avail], probs)
                bins <- mapply(
                    function(x, y) {
                        if(x == y) {
                            intersect(which(k == x), avail)
                        } else {
                            intersect(which(k >= x & k < y), avail)
                        }
                    },
                    quantiles[-length(quantiles)], quantiles[-1], SIMPLIFY = FALSE
                )
                
                VS <- replicate(ntrials, sapply(bins, sample, 1), simplify = FALSE) # stratified sample
            } else if(use_communities) {
                ## partition <- cluster_louvain(g)
                                        # try with a deterministic algorithm
                partition <- cluster_fast_greedy(g)

                sizes <- as.numeric(table(membership(partition)))
                pvec <- sqrt(sizes)/sum(sqrt(sizes))

                VS <- replicate(
                    ntrials,
                    as.numeric(select_by_comm_prob(n, g, partition, pvec, avail)),
                    simplify = FALSE
                )
            } else {
                VS <- replicate(ntrials, sample(avail, n), simplify = FALSE) # deg/knn constraints
            }
            
        } else {
            VS <- replicate(ntrials, sample(1:N, n), simplify = FALSE) # fully random
        }
    }

    KS <- lapply(VS, function(vs) k[vs])
    sorting <- lapply(KS, order)

    error <- mclapply(
        VS, function(x) obj_fn(x, y, Y, bparam, optimize_weights = optimize_weights),
        mc.cores = ncores, mc.silent = isFALSE(verbose)
    )

    dl <- vector("list", length(VS))
    for(i in seq_along(dl)) {
        dl[[i]] <- list(
            vs = VS[[i]][sorting[[i]]],
            error = error[[i]],
            ks = KS[[i]][sorting[[i]]]
        )
    }
    if(optimize_weights) {
        WS <- mclapply(VS, quadoptm, y = y, Y = Y, mc.cores = ncores, mc.silent = isFALSE(verbose))
        
        for(i in seq_along(dl)) dl[[i]]$ws <- WS[[i]][sorting[[i]]]
    }

    return(dl)
}

generate_data <- function(g, Y, y, n, ntrials = 100, optimize_weights = FALSE, verbose = FALSE) {
    random <- make_dataset(n, ntrials, doublewell_parms$Ds, y, Y)
    random_errors <- sapply(random, `[[`, "error")
    optimized <- make_dataset(
        n, ntrials, doublewell_parms$Ds, y, Y, optimize = TRUE, optimize_weights = optimize_weights,
        verbose = verbose
    )
    optimized_errors <- sapply(optimized, `[[`, "error")
    idx <- which.min(optimized_errors)
    ref <- optimized[[idx]]$vs
    fixingdeg <- make_dataset(n, ntrials, doublewell_parms$Ds, y, Y, comps = ref)
    fixingdeg_errors <- sapply(fixingdeg, `[[`, "error")

    rvs <- t(sapply(random, `[[`, "vs"))
    rks <- t(sapply(random, `[[`, "ks"))
    fvs <- t(sapply(fixingdeg, `[[`, "vs"))
    fks <- t(sapply(fixingdeg, `[[`, "ks"))
    ovs <- t(sapply(optimized, `[[`, "vs"))
    oks <- t(sapply(optimized, `[[`, "ks"))

    dl <- list(
        re = random_errors,
        fe = fixingdeg_errors,
        oe = optimized_errors,
        rvs = rvs,
        fvs = fvs,
        ovs = ovs,
        rks = rks,
        fks = fks,
        oks = oks
    )
    return(dl)
}

### Other useful functions
placelabel <- function(label, x, y, ...) { # defined in infile?
    xlim <- par("usr")[1:2]
    ylim <- par("usr")[3:4]
    xpos <- xlim[1] + x*(xlim[2] - xlim[1])
    ypos <- ylim[1] + y*(ylim[2] - ylim[1])
    if(par("xlog")) xpos <- 10^xpos
    if(par("ylog")) ypos <- 10^ypos
    text(xpos, ypos, label, ...)
}

#################
## Unused code ##
#################

## get_gcc <- function(g) {
##     if(is_directed(g)) {
##         comps <- components(g, mode = "weak")
##     } else comps <- components(g)

##     gcc_id <- which.max(comps$csize)
##     vids <- V(g)[comps$membership == gcc_id]
##     g <- induced_subgraph(g, vids)
##     return(g)
## }


## ### Subgraph analysis
## random_induced_subgraph <- function(vs, g) {
##     ki <- k[vs]
##     u <- sapply(ki, function(x) sample(as.numeric(V(g)[which(k == x)]), 1))
##     induced_subgraph(g, u)
## }

## n_random_subgraphs <- function(n, vs, g) {
##     lapply(seq(n), function(x) random_induced_subgraph(vs, g))
## }



## Methods above may not be appropriate for Kuramoto model
## Reason: a random node is near the mean---may not be a special node in any way.

## ### Kuramoto
## kuramoto <- function(t, θ, parms) {
##     with(parms, {# α, ω, θ, k
##         dθ <- ω - D*(colSums(A*t(outer(θ, θ, `-`)))/k)
##         return(list(c(dθ)))
##     })
## }

## kuramoto.state <- function(θ) abs(sum(exp(1i * θ))/N)

## solve_kuramoto <- function(θ, ω, Ds, A, times = seq(from = 0, to = 100, by = 1)) {
##     k <- rowSums(A)
##     results <- lapply(Ds, function(D) {
##         ode(θ, times, kuramoto, parms = list(D = D, ω = ω, k = k, A = A))
##     })

##     do.call(rbind, lapply(results, function(df) df[nrow(df), -1]))
## }

## kuramoto_parms <- list(
##     θ.init = function(N) 2*pi*runif(N),
##     ω = 1,# + runif(N, -.05, .01),
##     Ds = seq(0, .25, by = .0005),
##     αs = 1
## )

## The below is not well enough developed for use in this project.

## ### Wilson-Cowan Model of Neuron Excitation/Inhibition
## ## to implement the delay, need to solve with dede instead of ode and include a tau parameter
## wilsoncowan <- function(t, z, parms) { # x needs to be c(E, I)
##     with(parms, {
##         E <- z[1:N]
##         I <- z[(N+1):(2*N)]
        
##         dE <- -E + (S_Emax - E)*S_E(c1*E - c2*I + c5*colSums(A*E) + P)
##         dI <- -I + (S_Imax - I)*S_I(c3*E - c4*I)
##         return(list(c(dE, dI)))
##     })
## }

## solve_wilsoncowan <- function(E, I, c5s, A, N, times = 0:15) {
##     z <- c(E, I)
##     if("parallel" %in% sessionInfo()$basePkgs) {
##         ncores <- detectCores() - 1
##         results <- mclapply(
##             c5s, function(c5) {
##                 ode(z, times, wilsoncowan, parms = c(wilsoncowan_parms, c5 = c5, A = A, N = N))
##             }, mc.cores = ncores, mc.silent = TRUE
##         )
##     } else  {
##         results <- lapply(c5s, function(c5) {
##             ode(z, times, wilsoncowan, parms = c(wilsoncowan_parms, c5 = c5, A = A, N = N))
##         })
##     }

##     output <- do.call(rbind, lapply(results, function(df) df[nrow(df), -1]))
##     colnames(output) <- c(paste0("E", 1:N), paste0("I", 1:N))
##     return(output)
## }

## wilsoncowan_parms <- list(
##     E.init = 1,
##     I.init = 1,
##     c1 = 16,
##     c2 = 12,
##     c3 = 15,
##     c4 = 3,
##     ##c5 = .1,
##     P = 0,
##     ## τ = 8, # time constant, ignoring
##     ## τ_d = 10, # delay, ignoring
##     S_Emax = 1, # this is a guess based on the original paper
##     S_Imax = 1, # this is a guess based on the original paper
##     S_E = function(x, a_E = 1.3, θ_E = 4) {(1/(1 + exp(-a_E*(x - θ_E)))) - (1/(1 + exp(a_E*θ_E)))},
##     S_I = function(x, a_I = 2, θ_I = 3.7) {(1/(1 + exp(-a_I*(x - θ_I)))) - (1/(1 + exp(a_I*θ_I)))},
##     c5s = seq(0.001, 15, length.out = lout),
##     Ps = seq(0, 7, length.out = lout)# bifurcation parameter; stimulation to 1 node?; 1.25
## )


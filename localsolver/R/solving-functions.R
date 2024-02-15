#' @import stats parallel deSolve
library(stats)
library(parallel)
library(deSolve)

#' Functions to standardize simulations of ODE/SDEs
#'
#' Paired functions and lists of parameters to standardize using deSolve's ode() to compute ground truth simulations/solutions of several ODEs on various networks. 
#' @name ODEs
#' @param t An arbitrary time.
#' @param x The current state of the system.
#' @param params The list of model parameters.
#' @details Names without dots are functions which compute derivatives and are in deSolve's standard form. Each function can be used in a one-time way to compute the derivative of a system given an arbitrary time t, current state x, and parameters. More commonly, passed to deSolve's ode() as the model to be simulated.
#' Dot names are lists of standard model parameters. Pass the required coupling strength like `c(.name, list(D = D))`. The adjacency matrix is also required for solutions on networks; pass it in the same way. If analyzing a single variable, x should have length 1 and pass A = matrix(0).
#' @return A vector of derivatives
#' @examples
#' library(parallel)
#' ncores <- detectCores()-1
#' library(igraph)
#' library(deSolve)
#' g <- readRDS("../data/dolphin.rds")
#' N <- vcount(g)
#' A <- as_adj(g, "both", sparse = FALSE)
#' k <- degree(g)
#' times <- 0:15
#' params <- c(.doublewell, list(A = A))
#' control <- list(times = times, ncores = ncores)
#' X <- solve_in_range(params$Ds, "D", doublewell, rep(params$xinit.low, N), params, control, "ode")
NULL

#' @rdname ODEs
#' @export
.doublewell <- list(
    xinit.low = 1, xinit.high = 5,
    r = c(1, 3, 5),
    D = 0.05, Ds = seq(0, 1, length.out = lout),
    u = 0, us.up = seq(0, 5, length.out = lout), us.down = seq(0, -5, length.out = lout),
    sigma = 1e-2,
    lower.basin.limit = 2, upper.basin.limit = 4
)

#' @rdname ODEs
#' @export
doublewell <- function(t, x, params) {
    with(params, {
                                        # to keep the multiplication consistent and allow for directed
                                        # networks, form an N x N matrix in which each value of x is
                                        # repeated in each row. These are the x_j values for each x_i.
                                        # Use outer() to remain consistent across models.
                                        # Using outer() with all 1s, 1s first, makes a column of x_1,
                                        # a column of x_2, etc. Order matters for outer(). 
                                        # Use rowSums to sum across rows (the /i/s), aka over the
                                        # columns (the /j/s).
        coupling <- D*rowSums(A*outer(rep(1, length(x)), x))
        dx <- -(x - r[1])*(x - r[2])*(x - r[3]) + coupling + u
        return(list(c(dx)))
    })
}

#' @rdname ODEs
#' @export
.SIS <- list( # Only the up direction makes sense for this model (epidemic threshold)
    xinit.low = 0.01,
    mu = 1,
    D = 0.05, Ds = seq(0, 1, length.out = lout),
    u = 0,
    sigma = 1e-5,
    lower.basin.limit = 0.1 # epidemic threshold, name matches direction terminology for doublewell
)

#' @rdname ODEs
#' @export
SIS <- function(t, x, params) {
    with(params, {
        x[x < 0] <- 0 # zero is a floor here, and dx should only make it grow
        coupling <- D*rowSums(A*outer(1 - x, x))
        dx <- coupling - mu*x# + u # u is not physical for this model
        return(list(c(dx)))
    })
}

#' @rdname ODEs
#' @export
.genereg <- list( # Only the down direction makes sense for this model (going towards cell death)
    xinit.high = 2,
    B = 1, f = 1, h = 2,
    D = 1, Ds = seq(0, 1, length.out = lout),
    u = 0, us.down = seq(0, -0.5, length.out = lout),
    sigma = 1e-5,
    upper.basin.limit = 0.1 # marks cell death, name matches direction terminology for doublewell
)

#' @rdname ODEs
#' @export
genereg <- function(t, x, params) {
                                        # Currently allows for x < 0 if the bifurcation parameter is u
                                        # Set x = 0 if x < 0?
    with(params, {
        x[x < 0] <- 0
        coupling <- D*rowSums(A*outer(rep(1, length(x)), (x^h)/(1 + (x^h))))
        dx <- ifelse(
            x > 0,
            -B*(x^f) + coupling + u,
            0
        )
        return(list(c(dx)))
    })
}

#' @rdname ODEs
#' @export
.mutualistic <- list( # Up direction could be invasion, down collapse
    xinit.low = 0.001, xinit.high = 10,
    B = 0.1, K = 5, C = 1, Dtilde = 5, E = 0.9, H = 0.1,
    D = 0.05, Ds = seq(0, 3, length.out = lout),
    u = 0, us.up = seq(0, 0.5, length.out = lout), us.down = seq(0, -5, length.out = lout),
    sigma = 1e-3,
    lower.basin.limit = 1, upper.basin.limit = 4
)

#' @rdname ODEs
#' @export
mutualistic <- function(t, x, params) {
    with(params, {
        x[x < 0] <- 0
        coupling <- D*rowSums(A*(outer(x, x)/(Dtilde + outer(E*x, H*x, `+`))))
        dx <- ifelse(
            x > 0,
            B + x*(1 - (x/K))*((x/C) - 1) + coupling + u,
            0
        )
        return(list(c(dx)))
    })
}

#' Simulate ODEs
#' 
#' @description Simulate a system of ODEs across a range of a control parameter.
#'
#' @param rng The range of values the control parameter takes. 
#' @param varname The name of the control parameter. Must match a parameter in params list.
#' @param func The dynamical model, a function returning a derivative
#' @param initialvalue A vector or scalar of initial value(s).
#' @param params A list of function parameters, including the adjacency matrix if necessary.
#' @param control A list including some subset of nsamples, spacing, deltaT, times, ncores, or silent, depending on the application.
#' @param kind To solve with or without noise
#' @param allsamples Not currently used
#' 
#' @details Solves an ODE or SDE across a range of parameter values. If `kind = "ode"`, requires deSolve. 
#'
#' @return An L x N matrix of final values, where L is the length of `rng` and N is the length of `initialvalue`.
#' @examples
#' library(parallel)
#' ncores <- detectCores()-1
#' library(igraph)
#' library(deSolve)
#' g <- readRDS("../data/dolphin.rds")
#' N <- vcount(g)
#' A <- as_adj(g, "both", sparse = FALSE)
#' k <- degree(g)
#' times <- 0:15
#' params <- c(.doublewell, list(A = A))
#' control <- list(times = times, ncores = ncores)
#' X <- solve_in_range(params$Ds, "D", doublewell, rep(params$xinit.low, N), params, control, "ode")
#' @export
solve_in_range <- function(
                           rng, # the range of the parameter value, given as a vector
                           varname, # the name of the control parameter in the model, a string
                           func, # the dynamical model, returning a derivative with or without noise
                           initialvalue, # pass as vec or scalar... needs to be done outside
                           params = list(), # for the func
                           control = list(), # nsamples, spacing, deltaT, times, ncores, silent
                           kind = "ode", # or "sde"
                           allsamples = FALSE # returns 3D array
                           ) {
                                        # If mclapply won't run, try with ncores = 1
    if(!("ncores" %in% names(control))) control$ncores <- 1
    if(!("silent" %in% names(control))) control$silent <- FALSE # this doesn't seem to do much
    if(!("times" %in% names(control))) control$times <- seq(0, 10, by = 1)

    if(kind == "ode") {
        results <- mclapply(
            rng, function(val) {
                params[[varname]] <- val
                ode(initialvalue, control$times, func, params)
            }, mc.cores = control$ncores, mc.silent = control$silent
        )
                                        # Treat the final value of x_i as x_i^*
        result <- do.call(rbind, lapply(results, function(df) df[nrow(df), -1]))
        
    } else if(kind == "sde") {
        results <- mclapply(
            rng, function(val) {
                params[[varname]] <- val 
                sde(initialvalue, control$times, func, params, control)
            }, mc.cores = control$ncores, mc.silent = control$silent
        )
                                        # Treat the final value of x_i as x_i^*
        result <- do.call(rbind, lapply(results, function(df) df[nrow(df), -1]))
    }
    rownames(result) <- seq_len(nrow(result))
    return(result)
}

### SDEs

#' Determine the number of time steps for a simulation
#'
#' Uses the desired time span and the simulation time step to compute how many time steps are required for the simulation.
#' @param times a sequence of time steps, in user units
#' @param maxT The last time step (in user units)
#' @param minT The first time step (in user units)
#' @param deltaT The size of the timestep for simulation
#' @return A scalar
#' @examples
#' determine_ntimesteps(times = 0:10, deltaT = 0.01)
#' @export
determine_ntimesteps <- function(times = numeric(), maxT = 1, minT = 0, deltaT = 0.01) {
                                        # Requires times be at least length 2.
                                        # Not sure that's what I want.
    if(length(times > 1)) {
        minT <- times[1]
        maxT <- times[length(times)]
    }
    
    (maxT - minT)/deltaT
}
    
#' Preallocate noise for sde()
#'
#' Generates a ntimesteps x nnodes matrix of Gaussian random numbers to use for solving SDEs.
#'
#' @param sigma The standard deviation of the noise process
#' @param ntimesteps The number of simulation time steps (probably the output of determine_ntimesteps())
#' @param nnodes The number of SDEs (i.e., number of nodes in a network)
#' @return A matrix
#' @examples
#' ntimesteps <- determine_ntimesteps(times = 0:10, deltaT = 0.01)
#' W <- preallocate_noise(0.01, ntimesteps, 5)
#' @export
preallocate_noise <- function(sigma, ntimesteps, nnodes) {
                                        # hard-coding Gaussian noise
                                        # Would like to allow for Levy flights, eventually (power-law noise?)
                                        # See here: https://en.wikipedia.org/wiki/L%C3%A9vy_process
                                        # And here: https://en.wikipedia.org/wiki/L%C3%A9vy_flight
                                        # Looks like options are Wiener/Gaussian (below), Poisson, and Cauchy.
                                        # It looks like determining the matrix W would be more complicated for either a Poisson or Cauchy process, though potentially worth doing.
    if(length(sigma) == 1) {
        matrix(
            rnorm(nnodes*ntimesteps, mean = 0, sd = sigma), ncol = nnodes, nrow = ntimesteps
        )
    } else {
        stopifnot(length(sigma) == nnodes) # require σ_i ∀ i, but not necessarily unique
        do.call(
            cbind,
            lapply(sigma, function(s) rnorm(ntimesteps, mean = 0, sd = s))
        )
    }
}

#' Simulate stochastic differential equations on networks
#'
#' Simulate stochastic differential equations using the Euler-Maruyama method. Output mimics the output of deSolve's ode().
#' @param initialvalue The inital value of each variable/node, scalar or vector
#' @param times A sequence, at least the first and last time in user units, but e.g. 0:10 works and is easy.
#' @param func The dynamical model, in deSolve's format. Should return a deterministic derivative. 
#' @param parms A list. Using deSolve's naming convention, the model parameters (including the adjacency matrix, if using). Must include an element called "sigma" for the standard deviation of the noise process.
#' @param control A list. Must include an element called "deltaT".
#' @return A data frame
#' @details Returns a data frame that looks like deSolve's ode() output: a column of time steps, then a column of values at each timestep for each variable in the model.
#' @examples
#' library(igraph)
#' g <- sample_gnm(10, 20)
#' N <- vcount(g)
#' A <- as_adj(g, "both", sparse = FALSE)
#' X <- sde(rep(.doublewell$xinit.low, N), 0:10, doublewell, c(.doublewell, list(A = A)), list(deltaT = 0.01))
#' @export
sde <- function(initialvalue, times, func, parms = list(), control = list()) { # `parms` b/c deSolve
                                        # must have a noise strength
    stopifnot("sigma" %in% names(parms))
                                        # must have a Δt
    stopifnot("deltaT" %in% names(control))
                                        # must have adjacency matrix for coupled SDEs, but this doesn't enforce it
    if("A" %in% names(parms)) N <- nrow(parms$A) else N <- 1
    if("nstatevars" %in% names(parms)) nstatevars <- parms$nstatevars else nstatevars <- 1

    ntimesteps <- determine_ntimesteps(times, control$deltaT)
    W <- preallocate_noise(parms$sigma, ntimesteps, N*nstatevars)

    simtimes <- seq_len(ntimesteps)
    showtimes <- (simtimes - 1)*control$deltaT # Take a note of this: times arg is [minT, maxT). I think this makes sense: We start at 0, store our location at the initial point, then start moving. I want a trace since the beginning. Doesn't affect the output of solve_in_range(). 

    if(is.null(names(initialvalue))) varnames <- 1:N else varnames <- names(initialvalue)
    ##if(nstatevars > 1) {
        ##varnames <- 
    
    x <- initialvalue
    X <- matrix(0, nrow = ntimesteps, ncol = N*nstatevars)
    for(timestep in simtimes) {
        X[timestep, ] <- x
        x <- x +
                                        # This step doesn't use deSolve, only the deSolve format
                                        # i.e., functions which are written to be used with deSolve
                                        # the func() here is just a function that computes a next step
                                        # (though only tested with a derivative-computing function).
                                        # scaling for Δt for the derivative
            func(timestep, x, parms)[[1]]*control$deltaT + # from deSolve docs, derivative is [[1]]
                                        # scaling for Δt for the noise process
            W[timestep, ]*sqrt(control$deltaT)
    }

    df <- as.data.frame(cbind(showtimes, X))
    colnames(df) <- c("time", varnames)
    return(df)
}

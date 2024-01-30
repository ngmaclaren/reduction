if(getwd() != "/user/neilmacl/Documents/reduction/analysis") if(interactive()) setwd("./analysis")
library(parallel)
ncores <- detectCores()-1
library(igraph)
source("../src/functions.R")
                                        # will need to do this a lot
get_error <- function(dl) sapply(dl, `[[`, "error")

load("../data/celegans.rda")

dw.env <- new.env()
sis.env <- new.env()
mut.env <- new.env()
load("../data/sentinel-nodesets-celegans-dw.RData", dw.env)
load("../data/sentinel-nodesets-celegans-SIS.RData", sis.env)
load("../data/sentinel-nodesets-celegans-mutualistic.RData", mut.env)

with(dw.env, {
    source("../src/functions.R", local = TRUE)
    ntrials <- 100 # There are 5000 node sets. Down sample to 100 for now.
    opt <- opts[sample.int(n = 5000, size = ntrials)]
    bestopt <- opt[which.min(sapply(opt, `[[`, "error"))][[1]]
    fixed <- make_dataset(n, ntrials, bparam, y, Y, comps = bestopt$vs)
    rand <- make_dataset(n, ntrials, bparam, y, Y)
    constr <- make_dataset(n, ntrials, bparam, y, Y, use_connections = TRUE)
    quant <- make_dataset(n, ntrials, bparam, y, Y, use_connections = TRUE, use_quantiles = TRUE)
})
with(dw.env, lapply(list(opt, fixed, rand, constr, quant), function(dl) summary(get_error(dl))))
with(sis.env, {
    source("../src/functions.R", local = TRUE)
    ntrials <- 100 # There are 5000 node sets. Down sample to 100 for now.
    opt <- opts[sample.int(n = 5000, size = ntrials)]
    bestopt <- opt[which.min(sapply(opt, `[[`, "error"))][[1]]
    fixed <- make_dataset(n, ntrials, bparam, y, Y, comps = bestopt$vs)
    rand <- make_dataset(n, ntrials, bparam, y, Y)
    constr <- make_dataset(n, ntrials, bparam, y, Y, use_connections = TRUE)
    quant <- make_dataset(n, ntrials, bparam, y, Y, use_connections = TRUE, use_quantiles = TRUE)
})
with(sis.env, lapply(list(opt, fixed, rand, constr, quant), function(dl) summary(get_error(dl))))
with(mut.env, {
    source("../src/functions.R", local = TRUE)
    ntrials <- 100 # There are 5000 node sets. Down sample to 100 for now.
    opt <- opts[sample.int(n = 5000, size = ntrials)]
    bestopt <- opt[which.min(sapply(opt, `[[`, "error"))][[1]]
    fixed <- make_dataset(n, ntrials, bparam, y, Y, comps = bestopt$vs)
    rand <- make_dataset(n, ntrials, bparam, y, Y)
    constr <- make_dataset(n, ntrials, bparam, y, Y, use_connections = TRUE)
    quant <- make_dataset(n, ntrials, bparam, y, Y, use_connections = TRUE, use_quantiles = TRUE)
})
with(mut.env, lapply(list(opt, fixed, rand, constr, quant), function(dl) summary(get_error(dl))))

                                        # Get the reference node set info from the env and make it
                                        # available at top level
dl <- with(sis.env, list(opt = opt, fixed = fixed, rand = rand, constr = constr, quant = quant))
xdf <- data.frame(lapply(dl, get_error))
                                        # same, but for the comparisons
dw.y <- with(dw.env, y)
dw.Y <- with(dw.env, Y)
dw.bparam <- with(dw.env, bparam)
sis.y <- with(sis.env, y)
sis.Y <- with(sis.env, Y)
sis.bparam <- with(sis.env, bparam)
mut.y <- with(mut.env, y)
mut.Y <- with(mut.env, Y)
mut.bparam <- with(mut.env, bparam)

ydf <- data.frame(
    lapply(# go through the reference node sets
        dl,
        # for each, compare all node sets to the other dynamics
        function(l) sapply(l, function(l) obj_fn(l$vs, dw.y, dw.Y, dw.bparam))
    )
)

xlim <- range(unlist(xdf))
ylim <- range(unlist(ydf))

pdf("robustness-test-celegans-SIS_dw.pdf")
palette("Set 1")
plot(NULL, xlim = xlim, ylim = ylim, ylab = "Double-well", xlab = "SIS", main = "Error", log = "xy")
legend(
    "bottomright", col = 1:5, pch = 1, pt.lwd = 2, bty = "n",
    legend = c("Optimized", "Fixed-degree", "Random", "Constrained", "Quantiled")
)
for(i in seq_along(xdf)) points(xdf[[i]], ydf[[i]], col = adjustcolor(i, .5), pch = 1)
for(i in seq_along(xdf)) points(mean(xdf[[i]]), mean(ydf[[i]]), col = i, pch = 1, cex = 2, lwd = 3)
dev.off()

if(interactive()) {
    setwd("/user/neilmacl/Documents/reduction/analysis/")
    dynamics <- "genereg"
} else {
    library(optparse)
    optionlist <- list(
        make_option(
            c("-d", "--dynamics"), type = "character", default = "dw",
            help = "The dynamics to simulate on each network. Default is 'dw'. Options: 'dw', 'SIS', 'genereg', and 'mutualistic'."
        )
    )
    args <- parse_args(
        OptionParser(option_list = optionlist),
        convert_hyphens_to_underscores = TRUE
    )
    dynamics <- args$dynamics
}

                                        # Libraries and functions
library(parallel)
ncores <- detectCores()-1
library(igraph, lib.loc = "/user/neilmacl/rlocal/")
library(latex2exp, lib.loc = "/user/neilmacl/rlocal/")
library(sfsmisc)

                                        # Options
palette("Tableau 10")

infile <- paste0("../data/compare-methods-", dynamics, ".RData")
outfile <- paste0("../img/study-knn-", dynamics, ".pdf")
nets <- c(
    "dolphin", "celegans", "proximity", "euroroad", "email", "er", "ba", "hk", "gkk", "lfr"
)
envs <- paste(nets, "env", sep = "_")

load(infile)

with(email_env, {
    nbrmax <- function(v, g, k) {
        nbr <- neighbors(g, v)
        max(k[nbr])
    }
    nbrsd <- function(v, g, k) {
        nbr <- neighbors(g, v)
        sd(k[nbr])
    }
    vs <- do.call(rbind, lapply(opts, `[[`, "vs"))
    ks <- do.call(rbind, lapply(opts, `[[`, "ks"))
    knns <- matrix(knn[vs], ncol = n)
    lcc <- transitivity(g, "localundirected")
    lccs <- matrix(lcc[vs], ncol = n)
    nmx <- sapply(V(g), nbrmax, g, k)
    nmxs <- matrix(nmx[as.numeric(vs)], ncol = n)
    nsd <- sapply(V(g), nbrsd, g, k)
    nsds <- matrix(nsd[as.numeric(vs)], ncol = n)
    print(quantile(nmx, seq(0, 1, 0.05)))
    print(summary(as.numeric(nmxs)))
})


with(email_env, quantile(nmx, seq(0, 1, 1/20)))

with(email_env, {
    print(table(ks))
    print(summary(knns))
})

with(email_env, {
    print(summary(lcc))
    print(summary(as.numeric(lccs)))
})

with(email_env, print(quantile(k, probs = seq(0, 1, 1/20))))

with(ba_env, {
    print(summary(k))
    print(summary(as.numeric(ks)))
    print(quantile(k, probs = seq(0, 1, 1/20)))
})

with(email_env, {
    print(summary(knn))
    print(summary(as.numeric(knns)))
})

with(email_env, print(quantile(knn, probs = seq(0, 1, 1/40))))


with(email_env, print(head(ks)))
with(email_env, print(quantile(k, seq(0, 1, 1/n))))



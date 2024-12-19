library(sfsmisc)
library(optNS)
library(igraph)

networks <- c(
    "dolphin", "proximity", "celegans", "euroroad", "email", "gkk", "ba", "hk", "lfr", "er",
    "drosophila", "reactome", "route_views", "spanish", "foldoc", "tree_of_life", "word_assoc",
    "enron", "marker_cafe", "prosper"
)
Ns <- sapply(networks, function(net) vcount(readRDS(paste0("../data/", net, ".rds"))))

save_plots <- FALSE # TRUE

loc <- "../img/dynamics-comparisons/"

comps <- list() # this will be where we store the comparison results

for(network in networks) {
    ## network <- "tree_of_life"

    comps[[network]] <- list

    imgfile <- paste0(loc, network, ".pdf")

    Ylist <- readRDS(paste0("../data/fullstate-", network,".rds")) # still fixing dolphin network

    ns.types <- c("opt", "fixed", "rand")
    legendtext <- c("Optimized", "Degree-preserving", "Random")

    make_plotdata <- function(dynamicsA, dynamicsB, ns.types, Ylist) {
        A <- readRDS(paste0("../data/ns-", network, "_", dynamicsA, ".rds"))[ns.types]
        A.Y <- Ylist[[dynamicsA]]
        A.y <- rowMeans(A.Y)

        B <- readRDS(paste0("../data/ns-", network, "_", dynamicsB, ".rds"))[ns.types]
        B.Y <- Ylist[[dynamicsB]]
        B.y <- rowMeans(B.Y)

        list(
            x = lapply(A, get_error),
            y = lapply(A, function(ns) sapply(ns, function(x) obj_fn(x$vs, B.y, B.Y)))
        )
    }

    dyns <- c("doublewell", "mutualistic", "SIS", "genereg")
    conds <- expand.grid(dyns, dyns, stringsAsFactors = FALSE)
    colnames(conds) <- c("dynamicsB", "dynamicsA")
    labnames <- c("Double-well", "Mutualistic species", "SIS", "Gene-regulatory")
    labnames <- expand.grid(labnames, labnames, stringsAsFactors = FALSE)
    colnames(labnames) <- c("nameB", "nameA")
    conds <- cbind(conds, labnames)
    conds <- conds[conds$dynamicsA != conds$dynamicsB, ]
    conds <- conds[, c("dynamicsA", "dynamicsB", "nameA", "nameB")]

    comps[[network]] <- conds[, 1:2]
    comps[[network]]$count <- NA

    show_axis_labels <- TRUE # FALSE
    unit <- 2.5
    ht <- 4*unit
    wd <- 4*unit
    labelsize <- 1.25
    colors <- list(opt = "#3584e4", fixed = "#ff7800", rand = "#33d17a")
    pchs <- list(opt = 1, fixed = 2, rand = 0)

    if(save_plots) {
        pdf(imgfile, height = ht, width = wd)
    } else {
        dev.new(height = ht, width = wd)
    }

    par(mfcol = c(4, 4))
    if(show_axis_labels) {
        par(mar = c(4, 4.1, 1, 0.9), pty = "s")
    } else {
        par(mar = c(2.75, 2.75, 0.25, 0.25), pty = "s")
    }

    eqs <- list(
        eqno = 1:4,
        dynamics = dyns,
        name = c("Double-well", "Mutualistic species", "SIS", "Gene-regulatory"),
        model = c(
            expression(frac(d*x[i], d*t) == -(x[i] - r[1])*(x[i] - r[2])*(x[i] - r[3]) + D*sum(a[ij]*x, j=1,N)),
            expression(frac(d*x[i], d*t) == B + x[i]*bgroup("(", 1 - frac(x[i], K), ")")*bgroup("(", frac(x[i], C) - 1, ")") + D*sum(a[ij]*frac(x[i]*x[j], tilde(D) + E*x[i] + H*x[j]), j=1, N)),
            expression(frac(d*x[i], d*t) == -mu*x[i] + lambda*sum(a[ij]*(1 - x[i])*x[j], j=1, N)),
            expression(frac(d*x[i], d*t) == -B*x[i]^f + D*sum(a[ij]*frac(x[j]^h, 1 + x[j]^h), j=1, N))
        )
    )

    i <- 1
    eq <- 1
    for(pltno in seq(16)) {
        if(any(pltno == 1 + seq(0, 15, by = 5))) {
            plot(NULL, axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1))
            ## next
            text(0.5, 0.66, labels = eqs$name[eq], cex = 1.5, font = 2)
            text(0.5, 0.33, labels = eqs$model[eq], cex = switch(eq, 1, 0.75, 1, 1))
            eq <- eq + 1
            next
        }
        
        plotdata <- make_plotdata(conds[i, "dynamicsA"], conds[i, "dynamicsB"], ns.types, Ylist)
        x <- plotdata$x
        y <- plotdata$y
        
        xlim <- range(unlist(x))
        ylim <- range(unlist(y))
        comps[[network]][
            i, "count"
            ##comps$dynamicsA == conds[i, "dynamicsA"] & comps$dynamicsB == conds[i, "dynamicsB"], "count"
        ] <- mean(y$fixed) < mean(y$opt) # this is what I want to know.

        plot(NULL, xlim = xlim, ylim = ylim, xlab = "", ylab = "", xaxt = "n", yaxt = "n", log = "xy")
        eaxis(1, n.axp = 1, cex.axis = labelsize)
        eaxis(2, n.axp = 1, cex.axis = labelsize)

        if(show_axis_labels) {
            title(xlab = conds[i, "nameA"], cex.lab = labelsize)
            title(ylab = conds[i, "nameB"], cex.lab = labelsize)
        }
        for(nst in ns.types) points(x[[nst]], y[[nst]], col = adjustcolor(colors[[nst]], .25), pch = pchs[[nst]])
        for(nst in ns.types) points(mean(x[[nst]]), mean(y[[nst]]), col = colors[[nst]],
                                    pch = pchs[[nst]], lwd = 3, cex = 3)
        
        if(i == 1) {
            legend(
                "bottomright", bty = "n", col = unlist(colors), pch = unlist(pchs), pt.lwd = 2, pt.cex = 2,
                cex = labelsize, legend = legendtext
            )
        }

        i <- i + 1
    }

    if(save_plots) dev.off()

    ## library(latex2exp)
    ## ns <- readRDS(paste0("../data/ns-", network, "_doublewell.rds"))$opt # SIS
    ## S <- ns[[which.min(get_error(ns))]]
    ## g <- readRDS(paste0("../data/", network, ".rds"))
    ## print(V(g)$name[S$vs])


    ## Ds <- sdn::.doublewell$Ds
    ## ht2 <- 5
    ## wd2 <- 5
    ## ls2 <- 1.75

    ##                                     # on Double-well
    ## if(save_plots) {
    ##     pdf(paste0("../img/dynamics-comparison-doublewellbifplot-SI-", network, ".pdf"), height = ht2, width = wd2)
    ## } else {
    ##     dev.new(height = ht2, width = wd2)
    ## }
    ## par(mar = c(5, 5, 1, 1), pty = "s")
    ## ## matplot(
    ## ##     Ds, Ylist$doublewell, type = "l", lty = 1, lwd = 0.5, col = "#babdb6",
    ## ##     xlab = TeX(r"(D)", italic = TRUE), ylab = TeX(r"($x^*$)", italic = TRUE),
    ## ##     font.lab = 3, cex.lab = ls2, cex.axis = ls2
    ## ## )
    ## plot(
    ##     NULL, ylim = range(Ylist$doublewell[, S$vs]), xlim = range(Ds), xlab = "D", ylab = "x",
    ##     font.lab = 3, cex.lab = ls2, cex.axis = ls2
    ## )
    ## lines(Ds, rowMeans(Ylist$doublewell), lty = 1, lwd = 6, col = "#000000")
    ## matlines(Ds, Ylist$doublewell[, S$vs], lty = 1, lwd = 4, col = "#3584e4")
    ## lines(Ds, rowMeans(Ylist$doublewell[, S$vs]), lty = 1, lwd = 8, col = "#9141ac")
    ## if(save_plots) dev.off()

    ##                                     # on mutualistic
    ## if(save_plots) {
    ##     pdf(paste0("../img/dynamics-comparison-mutualisticbifplot-SI-", network, ".pdf"), height = ht2, width = wd2)
    ## } else {
    ##     dev.new(height = ht2, width = wd2)
    ## }
    ## par(mar = c(5, 5, 1, 1), pty = "s")
    ## ## matplot(
    ## ##     Ds, Ylist$mutualistic, type = "l", lty = 1, lwd = 0.5, col = "#babdb6",
    ## ##     xlab = TeX(r"(D)", italic = TRUE), ylab = TeX(r"($x^*$)", italic = TRUE),
    ## ##     font.lab = 3, cex.lab = ls2, cex.axis = ls2
    ## ## )
    ## with(list(Ds = seq(0, 3, length.out = 100)), {
    ##     plot(
    ##         NULL, ylim = range(Ylist$mutualistic[, S$vs]), xlim = range(Ds), xlab = "D", ylab = "x",
    ##         font.lab = 3, cex.lab = ls2, cex.axis = ls2
    ##     )
    ##     lines(Ds, rowMeans(Ylist$mutualistic), lty = 1, lwd = 6, col = "#000000")
    ##     matlines(Ds, Ylist$mutualistic[, S$vs], lty = 1, lwd = 4, col = "#3584e4")
    ##     lines(Ds, rowMeans(Ylist$mutualistic[, S$vs]), lty = 1, lwd = 8, col = "#9141ac")
    ## })
    ## if(save_plots) dev.off()

    ##                                     # on SIS
    ## if(save_plots) {
    ##     pdf(paste0("../img/dynamics-comparison-SISbifplot-SI-", network, ".pdf"), height = ht2, width = wd2)
    ## } else {
    ##     dev.new(height = ht2, width = wd2)
    ## }
    ## par(mar = c(5, 5, 1, 1), pty = "s")
    ## ## matplot(
    ## ##     Ds, Ylist$SIS, type = "l", lty = 1, lwd = 0.5, col = "#babdb6", ylim = c(0, 1),
    ## ##     xlab = TeX(r"(D)", italic = TRUE), ylab = TeX(r"($x^*$)", italic = TRUE),
    ## ##     font.lab = 3, cex.lab = ls2, cex.axis = ls2
    ## ## )
    ## plot(NULL, ylim = c(0, 1), xlim = range(Ds), xlab = "D", ylab = "x", font.lab = 3, cex.lab = ls2, cex.axis = ls2)
    ## lines(Ds, rowMeans(Ylist$SIS), lty = 1, lwd = 6, col = "#000000")
    ## matlines(Ds, Ylist$SIS[, S$vs], lty = 1, lwd = 4, col = "#3584e4")
    ## lines(Ds, rowMeans(Ylist$SIS[, S$vs]), lty = 1, lwd = 8, col = "#9141ac")
    ## if(save_plots) dev.off()

    ##                                     # on genereg
    ## if(save_plots) {
    ##     pdf(paste0("../img/dynamics-comparison-generegbifplot-SI-", network, ".pdf"), height = ht2, width = wd2)
    ## } else {
    ##     dev.new(height = ht2, width = wd2)
    ## }
    ## par(mar = c(5, 5, 1, 1), pty = "s")
    ## ## matplot(
    ## ##     Ds, Ylist$genereg, type = "l", lty = 1, lwd = 0.5, col = "#babdb6", ylim = c(0, 1),
    ## ##     xlab = TeX(r"(D)", italic = TRUE), ylab = TeX(r"($x^*$)", italic = TRUE),
    ## ##     font.lab = 3, cex.lab = ls2, cex.axis = ls2
    ## ## )
    ## plot(
    ##     NULL, ylim = range(Ylist$genereg[, S$vs]), xlim = range(Ds), xlab = "D", ylab = "x",
    ##     font.lab = 3, cex.lab = ls2, cex.axis = ls2
    ## )
    ## lines(Ds, rowMeans(Ylist$genereg), lty = 1, lwd = 6, col = "#000000")
    ## matlines(Ds, Ylist$genereg[, S$vs], lty = 1, lwd = 4, col = "#3584e4")
    ## lines(Ds, rowMeans(Ylist$genereg[, S$vs]), lty = 1, lwd = 8, col = "#9141ac")
    ## if(save_plots) dev.off()
}

for(i in seq_along(comps)) {
    comps[[i]]$network <- names(comps)[i]
    comps[[i]]$N <- Ns[i]
}

cdf <- do.call(rbind, comps)
rownames(cdf) <- seq(nrow(cdf))

with(list(agg = aggregate(count ~ network + N, data = cdf, FUN = sum, subset = network != "er")), {
    show(cor.test(log(agg$N), agg$count, method = "pearson"))
    ## cor.test(agg$N, agg$count, method = "spearman")
    ## cor.test(agg$N, agg$count, method = "kendall")
    if(save_plots) {
        pdf("../img/opt-fixed-comparison.pdf")
    } else {
        dev.new()
    }
    par(pty = "s", mar = c(5, 5, 1, 1))
    plot(
        count ~ N, data = agg, log = "x",
        pch = 1, cex = 2, lwd = 2, axes = FALSE,
        xlab = "", ylab = ""
    )
    box()
    eaxis(1, at = axTicks(1)[c(FALSE, TRUE)], cex.axis = 2)
    axis(2, cex.axis = 2)
    title(xlab = "N", cex.lab = 2, font.lab = 3)
    title(ylab = "Count", cex.lab = 2)
    if(save_plots) dev.off()
})

with(list(agg = aggregate(count ~ dynamicsA + dynamicsB, data = cdf, FUN = sum, subset = network != "er")), {
    dynamics <- c("doublewell", "mutualistic", "SIS", "genereg")
    dyns <- expand.grid(dynamics, dynamics, stringsAsFactors = FALSE)
    colnames(dyns) <- c("dynamicsA", "dynamicsB")
    ddf <- merge(dyns, agg, by = c("dynamicsA", "dynamicsB"), all.x = TRUE)
    mdf <- reshape(data = ddf, timevar = "dynamicsB", idvar = c("dynamicsA"), direction = "wide")
    rownames(mdf) <- mdf$dynamicsA
    mat <- as.matrix(mdf[, -1])
    colnames(mat) <- gsub("count.", "", colnames(mat))
    show(mat[dynamics, dynamics])
})

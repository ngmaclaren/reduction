library(optparse)
optionlist <- list(
    make_option(
        c("-s", "--save-plots"), action = "store_true", default = FALSE,
        help = "Store output figures as PDFs. If tabular output is created, saves to a text file. If FALSE [default], will display plot in an R graphics device (and silently fail to do so if graphics devices are not available.) Ignored for simulation files. Default is %default. "
    )
)

args <- parse_args(
    OptionParser(option_list = optionlist), args = c("-s"),
    convert_hyphens_to_underscores = TRUE
)

library(igraph)
library(latex2exp, lib.loc = "/user/neilmacl/rlocal/")
palette("Tableau 10")

saveplots <- args$save_plots
load("../data/knnfig-email-dw.RData")
save_plots <- saveplots

                                        # assign colors    
colorref <- data.frame(
    k = sort(unique(dat$k)),
    color = 1:length(unique(dat$k))
)

                                        # plot
if(save_plots) {
    pdf("../img/knn-fig-noweights.pdf", height = 7, width = 14)
} else {
    dev.new(height = 7, width = 14)
}
ticksize <- 1.75
labelsize <- 1.75
par(mar = c(4, 4, 1, 1)+0.5, mfrow = c(1, 2))
plot(
    jitter(dat$k, amount = 0.2), dat$error, pch = 1, cex = 1.5, lwd = 3,
    col = adjustcolor(colorref$color[match(dat$k, colorref$k)], alpha.f = 1),
    cex.axis = ticksize, cex.lab = labelsize, xlab = "Degree", ylab = "Approximation error"
)
points(as.numeric(k[ref$vs]), rep(ref$error, n),
       pch = 3, cex = 3, lwd = 3, col = "gray30")
vquants <- data.frame(
    v = as.numeric(V(g)),
    k = degree(g),
    bc = betweenness(g, directed = FALSE, normalized = TRUE),
    cc = closeness(g, mode = "all", normalized = TRUE),
    tr = transitivity(g, "local"),
    knn = knn(g)$knn
)
mtext("A", line = -2, adj = 0.02, cex = labelsize, font = 2)
datr <- aggregate(error ~ v, data = dat, FUN = mean)
df <- merge(vquants, datr, by = "v")
plot(
    error ~ knn, data = df, pch = 1, cex = 1.5, lwd = 3,
    col = adjustcolor(colorref$color[match(df$k, colorref$k)], alpha.f = 1),
    cex.axis = ticksize, cex.lab = labelsize,
    xlab = "Avgerage nearest neighbor degree", ylab = "Approximation error"
)
points(
    df$knn[df$v %in% ref$vs], df$error[df$v %in% ref$vs],
    col = "gray30", pch = 3, cex = 3, lwd = 3
)
mtext("B", line = -2, adj = 0.02, cex = labelsize, font = 2)
if(save_plots) dev.off()


                                        # weighted
load("../data/knnfig-email-dw-weighted.RData")
save_plots <- saveplots
colorref <- data.frame(
    k = sort(unique(dat$k)),
    color = 1:length(unique(dat$k))
)

if(save_plots) {
    pdf("../img/knn-fig-withweights.pdf", height = 7, width = 14)
} else {
    dev.new(height = 7, width = 14)
}
ticksize <- 1.75
labelsize <- 1.75
par(mar = c(4, 4, 1, 1)+0.5, mfrow = c(1, 2))
plot(
    jitter(dat$k, amount = 0.2), dat$error, pch = 1, cex = 1.5, lwd = 3,
    col = adjustcolor(colorref$color[match(dat$k, colorref$k)], alpha.f = 1),
    cex.axis = ticksize, cex.lab = labelsize, xlab = "Degree", ylab = "Error"
)
points(as.numeric(k[ref$vs]), rep(ref$error, n),
       pch = 3, cex = 3, lwd = 3, col = "gray30")
mtext("A", line = -2, adj = 0.02, cex = labelsize, font = 2)
vquants <- data.frame(
    v = as.numeric(V(g)),
    k = degree(g),
    bc = betweenness(g, directed = FALSE, normalized = TRUE),
    cc = closeness(g, mode = "all", normalized = TRUE),
    tr = transitivity(g, "local"),
    knn = knn(g)$knn
)
datr <- aggregate(error ~ v, data = dat, FUN = mean)
df <- merge(vquants, datr, by = "v")
plot(
    error ~ knn, data = df, pch = 1, cex = 1.5, lwd = 3,
    col = adjustcolor(colorref$color[match(df$k, colorref$k)], alpha.f = 1),
    cex.axis = ticksize, cex.lab = labelsize,
    xlab = "Avg. nearest neighbor degree", ylab = "Avg. error"
)
points(
    df$knn[df$v %in% ref$vs], df$error[df$v %in% ref$vs],
    col = "gray30", pch = 3, cex = 3, lwd = 3
)
mtext("B", line = -2, adj = 0.02, cex = labelsize, font = 2)
if(save_plots) dev.off()

library(nlme)
library(sfsmisc)

save_plots <- FALSE # TRUE

dl <- readRDS("../data/ns-fmri.rds")
for(i in seq_along(dl)) {
    dl[[i]]$testerror$pid <- rep(i, 100)
}

error <- data.frame(
    opt = sapply(dl, function(x) mean(x$testerror$opt)),
    rand = sapply(dl, function(x) mean(x$testerror$rand))
)
## error$id <- seq_len(nrow(error))

if(save_plots) {
    pdf("../img/fmri.pdf")
} else {
    dev.new()
}
par(mar = c(5.5, 5.5, 0.5, 0.5), pty = "s")
plot(NULL, xlab = "", ylab = "", xlim = range(error), ylim = range(error), axes = FALSE)
box()
eaxis(1, cex.axis = 1.5)
eaxis(2, cex.axis = 1.5)
title(xlab = "Random", cex.lab = 2)
title(ylab = "Optimized", cex.lab = 2, line = 4)
abline(a = 0, b = 1)
points(error$rand, error$opt, col = adjustcolor(1, .5))
if(save_plots) dev.off()

dfs <- lapply(dl, function(x) as.data.frame(x$testerror))
df <- do.call(rbind, dfs)

rdf <- reshape(
    df, varying = c("opt", "rand"), v.names = "error", timevar = "ns.type", times = c("opt", "rand"),
    direction = "long"
)
rownames(rdf) <- seq_len(nrow(rdf))
rdf$ns.type <- factor(rdf$ns.type, levels = c("rand", "opt"))

model <- lme(error ~ ns.type, random = ~ 1 | pid, data = rdf)
summary(model)
round(summary(model)$tTable, 3)

rdf$logerror <- log(rdf$error)
m2 <- lme(logerror ~ ns.type, random = ~ 1 | pid, data = rdf)
summary(m2)

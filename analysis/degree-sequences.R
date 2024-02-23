## Maybe log the y scale, but don't mess with log-ing the x scale

library(igraph)
library(optNS)
library(ggplot2)
library(ggthemes)
library(gridExtra)

prep_for_ggplot <- function(dl) {# pass which list of node sets, e.g. opts[[1]]
    x <- get_ks(dl)
    if(is.matrix(x)) {
        x <- as.data.frame(t(get_ks(dl)))
        colnames(x) <- paste0("node", seq_len(ncol(x)))
        y <- reshape(x, varying = colnames(x), v.names = "k", timevar = "node", direction = "long")
        rownames(y) <- seq_len(nrow(y))
        y$node <- factor(y$node, ordered = TRUE)
    } else {
        y <- data.frame(node = 1, k = x, id = seq_len(length(x)))
        y$node <- factor(y$node, ordered = TRUE)
    }

    return(y)
}

plotit <- function(df, net, orig = FALSE, tickeach = TRUE, logx = FALSE) {
    g <- get(net)
    binwidth <- switch(net, dolphin = 1, ba = 1) #0.1
    trans <- switch(net, dolphin = "identity", ba = "identity")#"log10")
    breaks <- switch(net, dolphin = seq_len(max(degree(g))), ba = waiver())

    if(orig) {
        plt <- ggplot(df, aes(k)) +
            geom_histogram(binwidth = binwidth, fill = "gray60", color = "black", linewidth = 0.25)
    } else {
        plt <- ggplot(df, aes(k, fill = node)) +
            geom_histogram(binwidth = binwidth, color = "black", linewidth = 0.25) +
            scale_fill_tableau(guide = "none")
    }
    ## if(tickeach) plt <- plt + scale_x_continuous(breaks = seq_len(max(degs)))
    ## if(logx) plt <- plt + scale_x_log10()
    plt +
        theme_classic() +
        theme(
            text = element_text(size = 15, family = "Nimbus Sans"),
            panel.background = element_rect(color = "black", linewidth = 1, fill = NA)
        ) +
        labs(x = "Degree", y = "Frequency") +
        coord_cartesian(xlim = range(degree(g))) +
        scale_x_continuous(breaks = breaks, trans = trans) +
        scale_y_continuous(trans = trans)
    
}
        
    
##     if(orig) {
##         ggplot(df, aes(k)) +
##             theme_classic() +
##             theme(
##                 text = element_text(size = 15, family = "Nimbus Sans"),
##                 panel.background = element_rect(color = "black", linewidth = 1, fill = NA)
##             ) +
##             geom_histogram(binwidth = 1, fill = "gray60")+ #, color = "black", linewidth = 0.25) +
##             labs(x = "Degree", y = "Frequency") +
##             coord_cartesian(xlim = c(1-0.25, max(degs)+0.25)) +
##             scale_x_continuous(breaks = seq_len(max(degs)))
##     } else {
##         ggplot(df, aes(k, fill = node)) +
##             theme_classic() +
##             theme(
##                 text = element_text(size = 15, family = "Nimbus Sans"),
##                 panel.background = element_rect(color = "black", linewidth = 1, fill = NA)
##             ) +
##             geom_histogram(binwidth = 1)+ #, color = "black", linewidth = 0.25) +
##             scale_fill_tableau(guide = "none") +
##             labs(x = "Degree", y = "Frequency") +
##             coord_cartesian(xlim = c(1-0.25, max(degs)+0.25)) +
##             scale_x_continuous(breaks = seq_len(max(degs)))
##     }
## }


dolphin <- readRDS("../data/dolphin.rds")
dolphin.N <- vcount(dolphin)
dolphin.k <- degree(dolphin)

ba <- readRDS("../data/ba.rds")
ba.N <- vcount(ba)
ba.k <- degree(ba)

dolphin.opts <- readRDS("../data/dolphin-1-to-12.rds")
ba.opts <- readRDS("../data/ba-4.rds")

dolphin.dfs <- lapply(dolphin.opts, prep_for_ggplot)
ba.dfs <- lapply(ba.opts, prep_for_ggplot)
dolphin.orig <- data.frame(k = dolphin.k)
ba.orig <- data.frame(k = ba.k)

plotlist <- c(
    list(plotit(dolphin.orig, "dolphin", orig = TRUE)),
    lapply(dolphin.dfs[1:4], plotit, net = "dolphin"),
    list(plotit(ba.orig, "ba", orig = TRUE)),
    lapply(ba.dfs, plotit, net = "ba")
)
              

dev.new(height = 10, width = 14)
do.call(grid.arrange, c(plotlist, ncol = 2, as.table = FALSE))

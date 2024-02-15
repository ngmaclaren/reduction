library(igraph)

setwd("../data/")

dolphin <- graph_from_data_frame(read.table("dolphin.txt"), directed = FALSE)
celegans <- graph_from_data_frame(read.table("celegans.txt"), directed = FALSE)
proximity <- graph_from_data_frame(read.table("infectious.txt"), directed = FALSE)
euroroad <- graph_from_data_frame(read.table("euro-road.txt"), directed = FALSE)
email <- graph_from_data_frame(read.table("arenas-email.txt"), directed = FALSE)

er <- graph_from_data_frame(read.table("er.txt"), directed = FALSE)
ba <- graph_from_data_frame(read.table("ba.txt"), directed = FALSE)
hk <- graph_from_data_frame(read.table("hk.txt"), directed = FALSE)
set.seed(12345)
gkk <- sample_fitness_pl(1000, 2500, 2.25)
lfr <- graph_from_data_frame(read.table("lfr.txt"), directed = FALSE)

networks <- c("dolphin", "celegans", "proximity", "euroroad", "email", "er", "ba", "hk", "gkk", "lfr")
filenames <- paste0(networks, ".rds")

for(i in seq_along(networks)) {
    g <- simplify(largest_component(get(networks[i])))
    saveRDS(g, filenames[i])
}

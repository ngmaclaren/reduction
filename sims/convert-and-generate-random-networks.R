library(igraph)

                                        # Erdős-Rényi graph
er <- graph_from_data_frame(read.table("../data/er.txt"), directed = FALSE)
## plot(er, vertex.size = 3, vertex.label = "")

                                        # Barabási-Albert graph
ba <- graph_from_data_frame(read.table("../data/ba.txt"), directed = FALSE)
## plot(ba, vertex.size = 3, vertex.label = "")

                                        # Holme-Kim graph
hk <- graph_from_data_frame(read.table("../data/hk.txt"), directed = FALSE)
## plot(hk, vertex.size = 3, vertex.label = "")

                                        # Lancichinetti-Fortunato-Radicchi graph
lfr <- graph_from_data_frame(read.table("../data/lfr.txt"), directed = FALSE)
## plot(lfr, vertex.size = 3, vertex.label = "")

                                        # Goh-Kahng-Kim node fitness graph with Cho et al correction
get_gcc <- function(g) { # In case pre-1.5 version of igraph...
    if(is_directed(g)) {
        comps <- components(g, mode = "weak")
    } else comps <- components(g)

    gcc_id <- which.max(comps$csize)
    vids <- V(g)[comps$membership == gcc_id]
    g <- induced_subgraph(g, vids)
    return(g)
}

set.seed(12345)
gkk <- get_gcc(sample_fitness_pl(1000, 2500, 2.25))
## plot(gkk, vertex.size = 3, vertex.label = "")

save(er, file = "../data/er.rda")
save(ba, file = "../data/ba.rda")
save(hk, file = "../data/hk.rda")
save(lfr, file = "../data/lfr.rda")
save(gkk, file = "../data/gkk.rda")

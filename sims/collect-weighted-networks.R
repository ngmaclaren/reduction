library(igraph)

setwd("../data/")

check <- function(g) {
    strengths <- strength(g, mode = "total", loops = FALSE, weights = E(g)$weight)
    cat("\n",
        "No. nodes:", vcount(g), "\n",
        "No. edges:", ecount(g), "\n",
        "Is directed:", is_directed(g), "\n",
        "Is connected:", is_connected(g), "\n",
        "Is simple:", is_simple(g), "\n",
        "Max. weight:", max(strengths), "\n",
        "Min. weight:", min(strengths), "\n"
    )
}

.plot <- function(g, widths = E(g)$weight) {
    plot(g, vertex.label = "", vertex.size = 5, edge.arrow.size = 0.5, edge.width = widths)
}

train_bombers <- graph_from_data_frame(read.csv("train-bombers.csv"), directed = FALSE)

foodweb_dry <- graph_from_data_frame(read.csv("foodweb-dry.csv"), directed = TRUE)
foodweb_dry <- as.undirected(foodweb_dry, mode = "collapse")

foodweb_wet <- graph_from_data_frame(read.csv("foodweb-wet.csv"), directed = TRUE)
foodweb_wet <- as.undirected(foodweb_wet, "collapse")

windsurfers <- graph_from_data_frame(read.csv("windsurfers.csv"), directed = FALSE)

highschool <- graph_from_data_frame(read.csv("highschool.csv"), directed = TRUE)
highschool <- as.undirected(highschool, "collapse")

macaques <- graph_from_data_frame(read.csv("macaques.csv"), directed = TRUE)
macaques <- as.undirected(macaques, "collapse")

residence_hall <- graph_from_data_frame(read.csv("residence-hall.csv"), directed = TRUE)
residence_hall <- as.undirected(residence_hall, "collapse")

flights <- graph_from_data_frame(read.csv("flights.csv"), directed = TRUE) # very skewed weight distribution
flights <- largest_component(as.undirected(flights, "collapse"))

celegans_neural <- graph_from_data_frame(read.csv("celegans-neural.csv"), directed = FALSE)
celegans_neural <- simplify(celegans_neural)

unicodelang <- graph_from_data_frame(read.csv("unicodelang.csv"), directed = FALSE)
unicodelang <- simplify(unicodelang)
unicodelang <- largest_component(unicodelang - E(unicodelang)[E(unicodelang)$weight == 0])

proximity_weighted <- graph_from_data_frame(read.csv("infectious.csv"), directed = FALSE)
proximity_weighted <- simplify(proximity_weighted)

networks <- c(
    "train_bombers", "foodweb_dry", "foodweb_wet", "windsurfers", "highschool", "macaques", "residence_hall",
    "flights", "celegans_neural", "unicodelang", "proximity_weighted"
)
filenames <- paste0(networks, ".rds")

for(i in seq_along(networks)) saveRDS(get(networks[i]), filenames[i])

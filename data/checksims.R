netnames <- c(
    "dolphin", "celegans", "proximity", "euroroad", "email", "er", "gkk", "ba", "hk", "lfr",
    "drosophila", "powergrid", "reactome", "route_views", "spanish", "foldoc", "tree_of_life",
    "word_assoc", "internet_as", "enron", "marker_cafe", "prosper"
)
filenames <- sapply(netnames, function(nomen) paste0("fullstate-", nomen, ".rds"))

fullstates <- lapply(filenames, readRDS)

lapply(fullstates, function(fs) sapply(fs, dim))

lapply(fullstates, function(fs) sapply(fs, range))

## need to do a nice pattern matching, regexp, to find all the right files
## sb "ns-*-first.rds|ns-*-last.rds" or something like that.
## I think list.files() accepts regexp patterns

netnames <- c("dolphin", "celegans", "proximity", "euroroad", "email", "er", "gkk", "ba", "hk", "lfr", "drosophila", "reactome", "route_views", "spanish", "foldoc", "tree_of_life", "word_assoc", "enron", "marker_cafe", "prosper")
dynnames <- c("doublewell", "SIS", "genereg", "mutualistic")

conds <- expand.grid(netnames, dynnames, stringsAsFactors = FALSE)

ns_names <- apply(conds, 1, paste, collapse = "_")

combine_ns_parts <- function(ns_name) {
    first_name <- paste0("ns-", ns_name, "-first.rds")
    last_name <- paste0("ns-", ns_name, "-last.rds")

    first_ns <- readRDS(first_name)
    last_ns <- readRDS(last_name)

    ns <- mapply(function(a, b) c(a, b), first_ns, last_ns, SIMPLIFY = FALSE)

    ns
}

for(ns_name in ns_names) {
    print(ns_name)
    filename <- paste0("ns-", ns_name, ".rds")
    saveRDS(combine_ns_parts(ns_name), filename)
}

## ## here's the example, which works. Just need to apply it to all sim conds
## t1 <- readRDS("ns-spanish_genereg-first.rds")
## t2 <- readRDS("ns-spanish_genereg-last.rds")
## newt <- mapply(function(a, b) c(a, b), t1, t2, SIMPLIFY = FALSE)

## check it
newt <- combine_ns_parts(ns_names[grep("prosper", ns_names)[3]])
class(newt)
length(newt)
lengths(newt)
lapply(newt, `[[`, 1)
lapply(newt, `[[`, 51)

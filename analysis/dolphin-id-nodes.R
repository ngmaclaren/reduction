## Remember to change GBB/DART objective functions
## Identify n ∈ {1, 2, 3, 4} for double-well and n = 4 for two other dynamics (a panel from the big fig)

library(igraph)

load("dolphin-demo.RData")

dw_nodes <- optNS::get_vs(solns)

## then here can just pick, but need the big plot first. So, need to reorganize it. 

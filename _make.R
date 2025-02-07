
source("_targets.R")
tar_make(as_job = TRUE)


tar_prune()
tar_glimpse()
tar_visnetwork()



library(ggraph)
library(tidygraph)

tar_network() %>%
      as_tbl_graph() %>%
      ggraph(layout = 'lgl') +
      geom_edge_link(aes(end_cap = label_rect(node2.name)),
                     arrow = arrow(length = unit(4, 'mm'),
                                   type = "closed", angle = 15)) +
      geom_node_label(aes(label = name, fill = type))

install.packages("igraph")
install.packages("ggraph")
install.packages("tidygraph")
install.packages("dplyr")

library(dplyr)

df <- data.frame(
  gene = c(
    "CYP1A1","CYP1A1","CYP1A1","CYP1A1",
    "CYP1B1","CYP1B1","CYP1B1","CYP1B1","CYP1B1",
    "NQO1","NQO1","NQO1","NQO1",
    "AKR1C1","AKR1C1","AKR1C1","AKR1C1"
  ),
  drug = c(
    "Omeprazole","Rifampicin","Quercetin","Benzo[a]pyrene",
    "Tamoxifen","Paclitaxel","Docetaxel","Resveratrol","Curcumin",
    "Mitomycin C","β-lapachone","Dicoumarol","Doxorubicin",
    "Progesterone","Indomethacin","Flufenamic acid","Doxorubicin"
  )
)

library(igraph)
library(ggraph)
library(tidygraph)

g <- graph_from_data_frame(df)

tg <- as_tbl_graph(g)

ggraph(tg, layout = "fr") +
  
  # edges
  geom_edge_link(alpha = 0.5, color = "gray50") +
  
  # nodes (genes vs drugs)
  geom_node_point(aes(color = name %in% df$gene), size = 5) +
  
  # labels
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  
  # styling
  theme_void() +
  ggtitle("Drug–Gene Interaction Network (TNBC Key Genes)") +
  
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "steelblue"),
                     labels = c("Drug", "Gene"),
                     name = "Node Type")
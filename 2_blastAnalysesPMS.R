# Library
library(tidyverse)
library(ggthemes)
library(cowplot)
options(width = 220)
setwd("D:/BIo/PMS_analysis/PMS_test4")
# Read data
df <- read_delim("allVOI_nucl.tab", delim = "\t", col_names = FALSE)
names(df) <- c("query", "subject", "identity", "alignment_length", "mismatches", "gaps", "q.start", "q.end", "s.start", "s.end", "evalue", "bit score")

# Transform data 
#' split column to have individual pieces of information
#' Done both for query and for subject
#' id is the threshold for identity
id <- 90
#' al_id is the threshold for alignment. The script will compare the length of the
#' contig and the lenght of the alignment. The alignment lenght should cover at least
#' the 90% of the contig lenght.
al_id <- 90
df1 <- df %>%
  mutate(qPool = str_replace(query, "(PMS*\\d+).*", "\\1")) %>%
  mutate(qTrinity = str_replace(query, ".*(TRINITY_\\w*_c\\d_g+\\d_i+\\d+|Contig\\d+).*", "\\1")) %>%
  mutate(qVirus = str_replace(query, ".*_i\\d+_(.*)_len=.*|.*Contig\\d+_(.*)_len=.*", "\\1\\2"),
         qVirus = str_replace_all(qVirus, "-", " ")) %>%
  mutate(qLength = str_replace(query, ".*len=(\\d*)", "\\1")) %>%
  mutate(sPool = str_replace(subject, "(PMS*\\d+).*", "\\1")) %>%
  mutate(sTrinity = str_replace(subject, ".*(TRINITY_\\w*_c\\d+_g\\d+_i\\d+|Contig\\d+).*", "\\1")) %>%
  mutate(sVirus = str_replace(subject, ".*_i\\d_(.*)_len=.*|.*Contig\\d+_(.*)_len=.*", "\\1\\2"),
         sVirus = str_replace_all(sVirus, "-", " ")) %>%
  mutate(sLength = str_replace(subject, ".*len=(\\d*)", "\\1")) %>%
  select(qPool:sLength, identity:`bit score`, query, subject) %>%
  mutate(self = case_when(query == subject ~ 1,
                          TRUE ~ 0),
         identOver = case_when(identity >= id ~ 1,
                               TRUE ~ 0),
         alignment_identity = case_when(alignment_length > 300 ~ 1,
                                        TRUE ~ 0)) %>%
  mutate(Group = NA)
###############################################################################
# Summarise data
###############################################################################
# Histogram of identity
p1 <- df1 %>%
  ggplot(aes(identity)) +
  geom_histogram(bins = 100) +
  theme_tufte() +
  scale_x_continuous(breaks=seq(0, max(df1$identity), 5)) +
  geom_vline(xintercept = id, linetype = "dashed", colour = "red")
p2 <- df1 %>%
  mutate(tmp = "A") %>%
  ggplot(aes(tmp, identity)) +
  geom_boxplot() +
  theme_void() +
  labs(title = "Identity distribution") +
  coord_flip()
plot_grid(p2, p1, ncol = 1)


###############################################################################
# Remove duplicates
###############################################################################
# Extract the full list of query IDs
qq <- df1 %>% 
  arrange(desc(as.numeric(qLength))) %>%
  pull(query) %>%
  unique()
df1_backup <- df1
# Remove the duplicates
for ( i in 1:nrow(df1)){
  tmp <- NULL
  tmp <- df1 %>%
    filter(query == qq[i]) %>%
    filter(alignment_identity != 0) %>%
    filter(identOver == 1) %>%
    filter(subject != qq[i]) %>%
    pull(subject)
  qq <- qq[!(qq %in% tmp)]
  df1 <- df1 %>%
    filter(!(query %in% tmp))
}

###############################################################################
# Filter data
###############################################################################
# Create groups
df2 <- df1 %>%
  group_by(query) %>%
  mutate(Group = ifelse(identity > id, query, NA)) %>%
  select(-(mismatches:evalue)) %>%
  ungroup

# Remove contigs not in groups and with alignment lenght below the threshold
#' I believe the alignment lenght below the threshold is some cases are still
#' there due to the fact that some of the subjects never became a query
df3 <- df2 %>%
  filter(!is.na(Group)) %>%
  filter(alignment_identity == 1) %>%
  ungroup

# Remove query duplicates based on alignment_length, identity and bit score
df4 <- df3 %>%
  group_by(Group) %>%
  top_n(1, alignment_length) %>%
  top_n(1, identity) %>%
  top_n(1, `bit score`) %>%
  ungroup

#extract list to txt univirus list --> Extract a fasta file from this list
univirus <- df4$query
univirus2 <- unique(univirus)
write(univirus2, "listavirusUni.txt")

###############################################################################
library(network)

# Node list
sources <- df3 %>%
  distinct(query) %>%
  rename(label = query)

subject <- df3 %>%
  distinct(subject) %>%
  rename(label = subject)

nodes <- full_join(sources, subject, by = "label")
nodes <- nodes %>% rowid_to_column("id")

# Edge list
per_route <- df3 %>%
  group_by(query, subject) %>%
  summarise(weight = mean(identity)) %>%
  ungroup()

edges <- per_route %>%
  left_join(nodes, by = c("query" = "label")) %>%
  rename(from = id)

edges <- edges %>%
  left_join(nodes, by = c("subject" = "label")) %>%
  rename(to = id)

edges <- select(edges, from, to, weight)

routes_network <- network(edges, vertex.attr = nodes, matrix.type = "edgelist", ignore.eval = FALSE)

plot(routes_network, vertex.cex = .5)

library(tidygraph)
library(ggraph)

routes_tidy <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)

ggraph(routes_tidy) + geom_edge_link() + geom_node_point() + theme_graph()

ggraph(routes_tidy, layout = "graphopt") +
  geom_node_point() +
  geom_edge_link(aes(width = weight), alpha = 0.8) +
  scale_edge_width(range = c(0.2, 2)) +
  theme_graph()
write_csv(df4, "df4.csv")

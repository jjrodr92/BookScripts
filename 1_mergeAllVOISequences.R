# Library
library(tidyverse)
library(DT)
library(fs)
setwd("your_workig_path") 

# Read data
filenames <- list.files(".", recursive = TRUE, pattern="SampleID*-Virus-reads\\.csv")
list.df <- map(filenames, read_csv, col_types = cols("TRINITY" = col_character()))

# Combine virus names
df <- list.df %>%
  purrr:::reduce(full_join) %>%
  mutate(Pool = rep(gsub("\\/.*", "", filenames), sapply(list.df, nrow))) %>%
  unite("Name", c("Pool", "Contig"), remove = FALSE) %>%
  mutate(Organism = gsub(" ", "-", Organism)) %>%
  unite("Final_name", c("Name", "Organism"), remove = FALSE) %>%
  mutate(Final_name = paste0(Final_name, "_len=", Length)) %>%
  arrange(Name)

#cahge 
df$Final_name <-gsub("-Virus-reads.csv", "_map2.fasta", df$Final_name)
df$Name <-gsub("-Virus-reads.csv", "_map2.fasta", df$Name)
df$Final_name <-gsub("SampleID-", "SampleID", df$Final_name)
df$Name <-gsub("SampleID-", "SampleID", df$Name)

# Save data
write_csv(df, "allVOI.csv")

# Read fasta sequences
library(Biostrings)
sequences <- list.files(".", recursive = TRUE, pattern="sampleID-significant2.fasta")
list.df2 <- map(sequences, readDNAStringSet)

df2 <- do.call(c, list.df2)

names(df2) <- paste(rep(gsub("\\/.*", "", sequences), sapply(list.df2, length)), gsub(" len=.*", "", names(df2)), sep = "_")

tmp <- tibble(tmp = names(df2))

tmp <- tmp %>% left_join(df %>% select(Final_name, Name), by = c("tmp" = "Name"))
names(df2) == tmp$tmp
names(df2) <- tmp$Final_name




#tmp <- paste(rep(gsub("\\/.*", "", sequences), sapply(list.df2, length)), gsub(" len=.*", "", names(df2)), sep = "_")
#tmp2 <- tibble(tmp)
#names(df2) <- tmp
#df2 <- df2[order(df2@ranges@NAMES),]
#names(df2) <- df$Final_name


writeXStringSet(df2, "allVOI.fasta")




# Library
library(tidyverse)
library(DT)
library(fs)
library(ggthemes)
library(testthat)
setwd("your_working_directory")
# Read data
filenames <- dir_ls(".", type = "file", regexp ="*sort\\.tab")
list.df <- map(filenames, read_delim, delim = "\t", skip = 1, col_names = FALSE)

# Combine virus reads
df <- do.call(rbind, list.df)     

# Create final table
#Create a table with colnames = Pools, row = contigs, values = number of reads
#Summarise all the reads together

#Summarise all the reads together
df1 <- df %>% 
  mutate(X5 = rep(str_extract(filenames, "SampleID*\\d+"),
                  sapply(list.df, nrow))) %>%
  select(Contig = X1,
         Reads = X3,
         Pool = X5) %>%
  spread(Pool, Reads) %>%
  mutate(AboveThr500 = rowSums(select(., `SampleID1`, `SampleID2`:`SampleID3`) > 500)) %>%
  arrange(desc(AboveThr500)) %>%
  mutate(Virus = str_replace(Contig, ".*_i\\d+_(.*)_len=.*|.*Contig\\d+_(.*)_len=.*", "\\1\\2"),
         Virus = str_replace_all(Virus, "-", " "),
         Length = str_replace(Contig, ".*len=(\\d*)", "\\1")) %>%
  select(Virus, Length, AboveThr500, `SampleID1`, `SampleID2`:`SampleID3`, Contig)

# Test
test_dir("tests/testthat/")

# Visualise the data
datatable(df1,
          rownames = FALSE,
          filter = "top",
          extensions = 'Buttons',
          options = list(dom = 'Bsftlip',
                         buttons = c('csv', 'copy'))
)

# Save table df1 
#' Contains only the  best contigs
write_csv(df1, "bestContigs_ReadMapped.csv")

# Plot data
df1 %>%
  gather(Pool, value, -Contig, -AboveThr500, -Virus, -Length) %>%
  mutate(Name = str_extract(Contig, "PMS*\\d+")) %>%
  ggplot(aes(Pool, value)) +
  geom_col(position = "dodge") +
  geom_col(aes(fill = Contig), position = "dodge") +
  # facet_wrap(~Pool) +
  coord_flip() +
  theme_tufte() +
  theme(legend.position = "none")

# # Most abundant virus
df1 %>%
  gather(Pool, value, -Contig, -AboveThr500) %>% filter(value == max(value))
# 
# # Merge this table with all blasts
# load("allVOI_blast.Rdata")
# 
df4 <- df1 %>%
  left_join(df3 %>%
              select(query, subject, self, identity, alignment_length, qPool:sVirus),
            by = c("Contig" = "query")) %>%
  unique

# # Save data
write_csv(df4, "allVOI_ReadMapped.csv")
save(df4, file = "allVOI_ReadMapped.RData")
# 
###############################################################################
# Add fasta sequence
library("Biostrings")
setwd("your_working_directory")
fastaFile <- readDNAStringSet("uniViruses1.fasta")
Contig = names(fastaFile)
Sequence = paste(fastaFile)
y <- data.frame(Contig, Sequence)

y$Contig <-gsub("sampleID-significant2.fasta", "", y$Contig)

df1$Contig <-gsub("sampleID-significant2.fasta", "", df1$Contig)

df1$Sequence <- y$Sequence[match(df1$Contig, y$Contig)]



# Write out

write_csv(df1, "bestContigs_ReadMapped.csv")



library("Biostrings")
setwd("D:/BIo/PMS_analysis/PMS_test4/tables")

df1 = read_csv ("hittablefinal.csv")

fastaFile <- readAAStringSet("uniViruses1Prot.fasta")
Contig = names(fastaFile)
Sequence = paste(fastaFile)
y <- data.frame(Contig, Sequence)

y$Contig2 <-gsub("-Virus-readsuni.csv", "", y$Contig)
y$Contig3 <-gsub("_len_", "_len=", y$Contig2)
df1$Contig <-gsub("-Virus-readsuni.csv", "", df1$qseqid)

df1$Sequence <- y$Sequence[match(df1$qseqid, y$Contig)]


fastaFile1 <- readAAStringSet("sequence.fasta")
Contig1 = names(fastaFile1)
Sequence1 = paste(fastaFile1)
z <- data.frame(Contig1, Sequence1)
z$contig3 <- str_split_fixed(z$Contig1, "=", 2)
z$contig4 <- trimws(x = z$contig3)
df1$AAHit3 <- z$Sequence[charmatch(df1$sseqid, z$contig4)]

# Write out

write_csv(df1, "hitstablefinal.csv")
write_csv(z, "Z.csv")

# descriptive statistics

mainly_alpha_seq <- read.csv('alpha_proteins_sequences.txt')

for (i in 1:nrow(mainly_alpha_seq)) {
  mainly_alpha_seq$length[i] <- nchar(as.character(mainly_alpha_seq$seq[i]))
}

hist(mainly_alpha_seq$length, breaks=20)
summary(mainly_alpha_seq$length)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 12.0   123.8   256.5   334.5   440.8  2166.0

intron_seq <- read.csv('../introns_final.txt', header=FALSE)
colnames(intron_seq) <- 'seq'

for (i in 1:nrow(intron_seq)) {
  intron_seq$length[i] <- nchar(as.character(intron_seq$seq[i]))
}

hist(intron_seq$length, breaks=100)
summary(intron_seq$length)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 11.0    323.8    987.0   3708.9   3405.2 187573.0
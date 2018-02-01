# load package to access Ensembl
library(biomaRt)
mart = useMart(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl')

### DATASET

# list of gene IDs
IDs <- read.csv('mart_export.csv')
# subset with gene type == 'protein_coding'
IDs <- IDs[IDs$Gene.type == 'protein_coding', ] # end up with 160,705 IDs
geneIDs <- IDs$Gene.stable.ID

# get individual exons: 523,351 observations
exons <- biomaRt::getSequence(id = geneIDs,
                              type = 'ensembl_gene_id',
                              seqType = 'gene_exon',
                              mart = mart)
exons$type <- 'exon'
colnames(exons) <- c('seq', 'geneID', 'type')
#write.csv(exons, 'data_exons.csv')
#read.csv('data_exons.csv')

# get un-separated sequences of introns and exons: 22,650 observations
exons_introns <- biomaRt::getSequence(id = geneIDs,
                                      type = 'ensembl_gene_id',
                                      seqType = 'gene_exon_intron',
                                      mart = mart) # 5' to 3'
colnames(exons_introns) <- c('exon_intron', 'geneID')
#write.csv(exons_introns, 'data_exons-introns.csv')
#read.csv('data_exons_introns.csv')

# match exons to exons_introns, replace exons with \n to be used as divider
for (exon in exons$exons) {
  for (i in 1:length(exons_introns$exon_intron)) {
    exons_introns[[1]][i] <- gsub(
      pattern = exon,
      replacement = ' ',
      x = exons_introns$exon_intron[i]
    )
  }
}

introns_vec <- as.vector(strsplit(exons_introns$exon_intron, ' '))
introns_df <- as.data.frame(as.matrix(unlist(introns_vec)))
introns_df$type <- 'intron'
colnames(introns_df) <- c('seq', 'type')
#write.csv(introns_df, 'data_introns.csv')

# combined dataset
exons_df <- subset(exons, select=-geneID)
data <- rbind(introns_df, exons_df)
data$gc <- 0
data$seq <- as.character(data$seq)

### FEATURE EXTRACTION

# GC content
for (i in 1:length(data$seq)) {
  ratio <- stringr::str_count(data$seq[i], 'G|C') / nchar(data$seq[i]) * 100
  data$gc[i] <- ratio
}
write.csv(data, 'data.csv')

mean(data$gc[data$type=='exon']) # 52.44306
mean(data$gc[data$type=='intron']) # 47.10053

# average disparity later?

### LINEAR REGRESSION
histogram(data$gc)

### CLASSIFICATION

set.seed(1)

# create 60/40 split
rows <- sample(nrow(data)) # shuffle row indices
data_shuffled <- data[rows, ] # randomly reorder rows
split <- round(nrow(data) * .6)
data_train <- data_shuffled[1:split, ]
data_test <- data_shuffled[(split+1):nrow(data), ]

# fit logistic regression model
model <- lm()
p <- predict(model, data_test)
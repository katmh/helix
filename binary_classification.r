# load package to access Ensembl
library(biomaRt)

mart = useMart(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl')

# list of gene IDs
IDs <- read.csv('mart_export.csv')
# subset with gene type == 'protein_coding'
IDs <- IDs[IDs$Gene.type == 'protein_coding', ] # end up with 160,705 IDs

geneIDs <- IDs$Gene.stable.ID

# get individual exons
exons <- biomaRt::getSequence(id = geneIDs,
                              type = 'ensembl_gene_id',
                              seqType = 'gene_exon',
                              mart = mart)
colnames(exons) <- c('exon', 'geneID')
write.csv(exons, 'data_exons.csv')

# get un-separated sequences of introns and exons
exons_introns <- biomaRt::getSequence(id = geneIDs,
                                      type = 'ensembl_gene_id',
                                      seqType = 'gene_exon_intron',
                                      mart = mart) # 5' to 3'
colnames(exons_introns) <- c('exon_intron', 'geneID')
write.csv(exons, 'data_exons_introns.csv')

# match exons to exons_introns and replace with \n to be used as divider
sep_introns <- data.frame('', c('', ''))
for (exon in exons$exons) {
  for (i in 1:length(exons_introns$exon_intron)) {
    sep_introns[[i]] <- gsub(
      pattern = exon,
      replacement = '\n',
      x = exons_introns$exon_intron[i]
    )
  }
}

seq <- gsub(
  pattern = exon,
  replacement = '\n',
  x = seq
)

exons_introns$gene_exon_intron <- gsub(pattern = exons$gene_exon[1], replacement = '\n', x = exons_introns$gene_exon_intron)

introns <- as.vector(strsplit(exons_introns$gene_exon_intron, '\n'))

fake_ie <- c('abcdefg', 'abcgfedc')
fake_e <- c('abc', 'fe')

(fake_ie[1] <- gsub(
  pattern = fake_e[1],
  replacement = ' ',
  x = fake_ie[1]
))

# beware of variable scope
new_seq <- data.frame('')
for (exon in fake_e) {
  for (i in 1:length(fake_ie)) {
    (new_seq[[i]] <- gsub(
      pattern = exon,
      replacement = ' ',
      x = fake_ie[i]
    ))
  }
}

new_seq

fake_ie

str1 = 'text\nmoretext'
strsplit(str1, '\n')

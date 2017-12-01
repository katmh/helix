setwd('~/helix')

combinedAlphas <- read.csv('combinedAlphas.csv', header=TRUE, sep=',')

### separate sequence and secondary structure

# create dataset of separate seq and SS
ssSeqSep <- data.frame(do.call('rbind', strsplit(as.character(combinedAlphas$seqSS), '#', fixed=TRUE)))

# remove original seqSS
combinedAlphas <- subset(combinedAlphas, select=-seqSS)

# add separate dataset to original
combinedAlphas <- merge(combinedAlphas, ssSeqSep, by='row.names')
# remove extra added column; TODO: fix this later
combinedAlphas <- subset(combinedAlphas, select=-Row.names)

maxdisparities <- readLines('maxdisparities.txt')

hist(as.numeric(maxdisparities),
     breaks = 47,
     xlim = c(0,2.5))
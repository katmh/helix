# Set working directory
setwd('~/helix')

# Load data
SSversusDisp <- read.csv('SSversusDisp.txt', header=TRUE, sep=',')
# 3,835,872 residues

# Load ggplot2 package
library(ggplot2)

# Subset of residues classified as alpha helix (H), beta strand/ladder (E), or unclassified
# TODO: what is beta bridge (B)?
HEU <- subset(SSversusDisp, ss=='H' | ss=='E' | ss==' ')
# 2,931,644 residues

# Save as PNG
png('SSversusDisp-HEU.png',
    width=1000,
    height=700,
    units='px',
    bg='transparent')

ggplot(HEU,
       aes(x=disp, fill=ss)) +
  geom_density(alpha=.3) +
  xlab('Disparity') +
  ylab('Density') +
  scale_fill_discrete(name='Secondary Structure',
                      labels=c('Alpha Helix', 'Beta Strand', 'Unclassified'))
dev.off()

# Both of these plot ss type vs count
plot(x = SSversusDisp$ss)
ggplot(data = SSversusDisp,
       aes(x=ss, fill=ss)) +
  geom_bar()

maxDisp <- read.csv('max-disparities.txt', header = TRUE)

hist(x = maxDisp$disp,
     breaks = 78,
     xlab = 'Maximum Disparity',
     main = 'Distribution of Maximum Disparity Values')

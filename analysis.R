# Set working directory
setwd('~/helix')

# Load data
SSversusDisp <- read.csv('SSversusDisp.txt', header=TRUE, sep=',')

# Load ggplot2 package
library(ggplot2)

ggplot(SSversusDisp, aes(x=disp, fill=ss)) +
  geom_density(alpha=.3)
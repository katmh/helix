# set working directory
setwd('~/helix2/helix')

# normality of mean/max disparity of sequences sorted by sec struc
# conclusion: none of samples are normally distributed

a <- read.csv('mm_mainly_alpha.csv')
hist(a$mean) # possibly normal by inspection
hist(a$max) # very skewed left
# test normality
shapiro.test(a$mean[1:5000]) # p<2.2e-16; null hypothesis that it's normal is rejected
shapiro.test(a$max[1:5000]) # also p<2.2e-16

b <- read.csv('mm_mainly_beta.csv')
hist(b$mean) # possibly normal by inspection
hist(b$max) # significantly skewed left
# test normality
shapiro.test(b$mean[1:5000]) # p<2.2e-16; null hypothesis that it's normal is rejected
shapiro.test(b$max[1:5000]) # p<1.85e-11

ab <- read.csv('mm_alpha_beta.csv')
hist(ab$mean) # possibly normal by inspection
hist(ab$max) # significantly skewed left
# test normality
shapiro.test(ab$mean[1:5000]) # p<2.2e-16; null hypothesis that it's normal is rejected
shapiro.test(ab$max[1:5000]) # p<2.2e-16
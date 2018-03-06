library(dplyr)
library(caret)
library(ROCR)
library(plotROC)
library(ggplot2)
library(plotly)

### CLEAN

a <- read.csv('alpha_coding_filtered.csv', header=FALSE)
b <- read.csv('beta_coding_filtered.csv', header=FALSE)
ab <- read.csv('ab_coding_filtered.csv', header=FALSE)
fSS <- read.csv('fewSS_coding_filtered.csv', header=FALSE)

a$ss <- 'a'
b$ss <- 'b'
ab$ss <- 'ab'
fSS$ss <- 'fSS'

data <- rbind(a,b,ab,fSS)
data$code <- 1
colnames(data) <- c('seq', 'ss', 'code')
write.csv(data, 'data.csv')

# quick start
setwd('C:/Users/Jan_Huang/Desktop/helix_grindtime/by_secondary_structure')
data <- read.csv('data.csv')
data <- subset(data, select=c('seq', 'length', 'code', 'disp', 'GC', 'ss', 'maxORF'))
data$code <- as.factor(data$code)

# get sequences for beta, ab, and fewSS
b <- read.csv('beta_coding.csv')
ab <- read.csv('ab_coding.csv')
fSS <- read.csv('fewSS_coding.csv')

# make results reproducible
set.seed(500)

# load datasets
coding <- read.csv('coding_final.csv')
noncoding <- read.csv('noncoding_final.csv')

# compare sequence lengths
for (i in 1:nrow(coding)) {
  coding$length[i] <- nchar(as.character(coding$seq[i]))
}
for (i in 1:nrow(noncoding)) {
  noncoding$length[i] <- nchar(as.character(noncoding$seq[i]))
}

hist(coding$length, main='Length of Coding Sequences', xlim=c(0, 100000), breaks=100, xlab='Length')
hist(noncoding$length, main='Length of Noncoding Sequences', xlim=c(0, 100000), breaks=100, xlab='Length')

# remove sequences <45 bp or >50000 bp
coding <- coding[coding$length < 50000, ] # removes 1 outlier (length>100000)
noncoding <- noncoding[noncoding$length >= 45, ] # removes 9 sequences

# make datasets the same size by randomly choosing sequences from bigger dataset
fraction <- nrow(coding) / nrow(noncoding)
nc_index <- sample(x = 1:nrow(noncoding), size = round(fraction * nrow(noncoding)))
noncoding <- noncoding[nc_index, ]

# combine datasets and save as CSV for Python manipulation
coding <- subset(coding, select=c('seq', 'length'))
coding$code <- 1
noncoding$code <- 0
data <- rbind(coding, noncoding)
write.csv(data, 'data.csv')

### FEATURE EXTRACTION

# in Python I'll calculate GC content, disparity
# then I'll import the values here
disp <- as.matrix(read.csv('data_disp.csv'))
for (i in 1:nrow(disp)) {
  data$disp[i] <- disp[i]
}
GC <- as.matrix(read.csv('data_GC.csv'))
for (i in 1:nrow(GC)) {
  data$GC[i] <- GC[i]
}
write.csv(data, 'data.csv')

b_disp <- as.matrix(read.csv('beta_coding_disp.csv'))
for (i in 1:nrow(b_disp)) {
  b$disp[i] <- b_disp[i]
}
b_GC <- as.matrix(read.csv('beta_coding_GC.csv'))
for (i in 1:nrow(b_GC)) {
  b$GC[i] <- b_GC[i]
}

ab_disp <- as.matrix(read.csv('ab_coding_disp.csv'))
for (i in 1:nrow(ab_disp)) {
  ab$disp[i] <- ab_disp[i]
}
ab_GC <- as.matrix(read.csv('ab_coding_GC.csv'))
for (i in 1:nrow(ab_GC)) {
  ab$GC[i] <- ab_GC[i]
}

fSS_disp <- as.matrix(read.csv('fewSS_coding_disp.csv'))
for (i in 1:nrow(fSS_disp)) {
  fSS$disp[i] <- fSS_disp[i]
}
fSS_GC <- as.matrix(read.csv('fewSS_coding_GC.csv'))
for (i in 1:nrow(fSS_GC)) {
  fSS$GC[i] <- fSS_GC[i]
}

b$ss <- 'b'
ab$ss <- 'ab'
fSS$ss <- 'fSS'
data$ss <- ''
data[data$code==1,]$ss <- 'a'

data_other_ss <- rbind(b,ab,fSS)
for (i in 1:nrow(data_other_ss)) {
  data_other_ss$length[i] <- nchar(as.character(data_other_ss$seq[i]))
}
data_other_ss$code <- 1
write.csv(data_other_ss, 'data_other_ss.csv')

data <- rbind(data, data_other_ss)

ORF_lengths <- as.matrix(read.csv('data_ORF.csv'))
for (i in 1:nrow(data)) {
  data$maxORF[i] <- ORF_lengths[i]
}

### DESCRIPTIVE STATISTICS FOR FEATURES 

mean(data[data$code == '1', ]$GC) # 51.86%
mean(data[data$code == '0', ]$GC) # 46.49%
sd(data[data$code == '1', ]$GC) # 7.93%
sd(data[data$code == '0', ]$GC) # 12.02%

mean(data[data$code == '1', ]$disp) # .73
mean(data[data$code == '0', ]$disp) # .61
sd(data[data$code == '1', ]$disp) # .11
sd(data[data$code == '0', ]$disp) # .13

mean(data[data$code == '1', ]$maxORF) # 439.11
mean(data[data$code == '0', ]$maxORF) # 221.08
sd(data[data$code == '1', ]$maxORF) # 394.13
sd(data[data$code == '0', ]$maxORF) # 230.42

shapiro.test(data[data$code=='1',]$GC) # not normal; p=1.028e-07
shapiro.test(data[data$code=='0',]$GC) # not normal; p=6.503e-05
shapiro.test(data[data$code=='1',]$disp) # not normal; p=2.241e-07
shapiro.test(data[data$code=='0',]$disp) # not normal; p=2.2e-16
# conclusion: none are normal

fligner.test(data[data$code=='1',]$GC, data[data$code=='0',]$GC) # p=.405, variances are same
fligner.test(data[data$code=='1',]$disp, data[data$code=='0',]$disp) # p=.492, variances are same
# conclusion: variances are same

# Wilcoxon Rank Sum Test (compare mean of 2 samples)
# alternative: true location shift is not equal to 0
wilcox.test(data[data$code=='1',]$GC, data[data$code=='0',]$GC) # p=9.46e-14
wilcox.test(data[data$code=='1',]$disp, data[data$code=='0',]$disp) # p=2.2e-16
# conclusion: from diff distributions

ggplot(data=data, aes(data$disp, fill=code)) + geom_histogram(binwidth = .05) + labs(title='Histogram of Median Disparity', x='Median Disparity', y='Count') + scale_fill_discrete(name = '', labels = c('Noncoding', 'Coding'))

ggplot(data=data, aes(data$GC, fill=code)) + geom_histogram(binwidth = 1.15) + labs(title='Histogram of GC Percentage', x='GC Percentage', y='Count') + scale_fill_discrete(name = '', labels = c('Noncoding', 'Coding'))

ggplot(data=data, aes(data$maxORF, fill=code)) + geom_histogram() + labs(title='Histogram of Max ORF Length', x='Max ORF Length', y='Count') + scale_fill_discrete(name = '', labels = c('Noncoding', 'Coding'))

# visualize correlations btwn/among features
ggplot(data=data, aes(x=data$GC, y=data$disp)) + geom_point(aes(GC, disp, color=code)) + labs(title='Scatterplot of GC Percentage vs Median Disparity')

ggplot(data=data, aes(x=data$maxORF, y=data$disp)) + geom_point(aes(GC, disp, color=code)) + labs(title='Scatterplot')

### LOGISTIC REGRESSION

# 75/25 train/test split
index <- sample(x = 1:nrow(data), size = round(0.75 * nrow(data))) # indices of observations to be included in train set, drawing w/o replacement by default
train <- data[index,] # 1486 observations
test <- data[-index,] # rows that train didn't include; 495 observations

# 10-fold cross validation

# train model on GC + disp
model <- glm(code ~ GC + disp, family = 'binomial', train)
p <- predict(model, test, type='response')
summary(p)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0003895 0.3029141 0.4716623 0.4928045 0.6545059 0.9922203
p_class <- ifelse(p > .50, 1, 0)
confusionMatrix(p_class, test[['code']])
# accuracy: 75.08%
# 95% CI: (69.9%, 79.8%)
# sens: 76.51%
# spec: %73.47

# train model on GC, disp, and maxORF
model <- glm(code ~ GC + disp + maxORF, family = 'binomial', train)
p <- predict(model, test, type='response')
summary(p)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.006836 0.899453 0.954246 0.915833 0.982934 1.000000
p_class <- ifelse(p > .50, 1, 0)
confusionMatrix(p_class, test[['code']])
# accuracy: 75.08%
# 95% CI: (69.9%, 79.8%)
# sens: 76.51%
# spec: %73.47

### ROC

pred <- prediction(p, test$code)
sens <- performance(pred, measure='sens', x.measure='cutoff') # true positive rate
spec <- performance(pred, measure='spec', x.measure='cutoff') # true negative rate
plot(sens, main='Sensitivity and Specificity vs Cutoff', col='red', ylab='Sensitivity/Specificity')
plot(spec, add=TRUE, col='blue')
# work out geom_roc later

### FEATURE SELECTION

# remove redundant features
(correlationMatrix <- cor(data[,4:5]))
(highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.5))
# disp and GC have .125 correlation

# rank featurs by importance, learning vector quantization (LVQ)
importance <- varImp(model, scale=FALSE)
print(importance)
#       Overall
#GC    5.570123
#disp 12.255953

# recursive feature elimination (RFE), which uses random forests cross-validation
# this part is a work in progress, haven't fully interpreted it yet
control <- rfeControl(functions=rfFuncs, method='cv', number=10)
results <- rfe(data[,3:7], data[,2], sizes=c(3:7), rfeControl=control) # got 31 warnings about <=5 unique values
print(results)
predictors(results)
plot(results, type=c('g', 'o'))

### RANDOM FORESTS
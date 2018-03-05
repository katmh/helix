library(dplyr)
library(caret)
library(ROCR)
library(plotROC)
library(ggplot2)

# quick start
setwd('C:/Users/Jan_Huang/Desktop/helix_grindtime/by_secondary_structure')
data <- read.csv('data.csv')
data <- subset(data, select=c('seq', 'length', 'code', 'disp', 'GC'))
data$code <- as.factor(data$code)

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

### DESCRIPTIVE STATISTICS FOR FEATURES 

ggplot(data=data, aes(data$disp, fill=code)) + geom_histogram(binwidth = .05) + labs(title='Histogram of Median Disparity', x='Median Disparity', y='Count') + scale_fill_discrete(name = '', labels = c('Noncoding', 'Coding'))

ggplot(data=data, aes(data$GC, fill=code)) + geom_histogram(binwidth = 1.15) + labs(title='Histogram of GC Percentage', x='GC Percentage', y='Count') + scale_fill_discrete(name = '', labels = c('Noncoding', 'Coding'))

# visualize correlations btwn/among features
ggplot(data=data, aes(x=data$GC, y=data$disp)) + geom_point(aes(GC, disp, color=code)) + labs(title='Scatterplot of GC Percentage vs Median Disparity')

### LOGISTIC REGRESSION

# 75/25 train/test split
index <- sample(x = 1:nrow(data), size = round(0.75 * nrow(data))) # indices of observations to be included in train set, drawing w/o replacement by default
train <- data[index,] # 1486 observations
test <- data[-index,] # rows that train didn't include; 495 observations

# 10-fold cross validation

# train model
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
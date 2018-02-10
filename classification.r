library(caret)
library(ROCR)

set.seed(500)

setwd('C:/Users/Jan_Huang/Desktop/helix_grindtime')
data <- read.csv('dataset.csv')

### SIMPLE FEATURES ONLY

# 75/25 train/test split
index <- sample(x = 1:nrow(data), size = round(0.75 * nrow(data))) # indices of observations to be included in train set, drawing w/o replacement by default
train <- data[index,]
test <- data[-index,] # rows that train didn't include

# logistic regression
model <- glm(ie ~ gc + nstart + maxnstop, family = 'binomial', train)
p <- predict(model, test, type='response')
summary(p)
p_class <- ifelse(p > .50, 1, 0)
confusionMatrix(p_class, test[['ie']])
pred.A <- prediction(p, test$ie)
ROC.A <- performance(pred.A, measure='tpr', x.measure='fpr')

### SIMPLE FEATURES + DISPARITY (EXPERIMENTAL)

# import maxdisp
data <- read.csv('dataset.csv')

maxdisps <- as.matrix(read.delim('maxdisparities.txt', header=FALSE))
data$maxdisp <- maxdisps
data <- data[!(is.na(data$maxdisp)), ] # remove rows with NA in maxdisp column
#write.csv(data, 'dataset_disp.csv')

# 75/25 train/test split
index <- sample(x = 1:nrow(data), size = round(0.75 * nrow(data))) # indices of observations to be included in train set, drawing w/o replacement by default
train <- data[index,]
test <- data[-index,] # rows that train didn't include

# logistic regression w/ maxdisps
model <- glm(ie ~ gc + nstart + maxnstop + maxdisp, family = 'binomial', train)
p <- predict(model, test, type='response')
summary(p)
p_class <- ifelse(p > .50, 1, 0)
table(p_class)
table(p_class, test[['ie']])
confusionMatrix(p_class, test[['ie']])

pred.B <- prediction(p, test$ie)
ROC.B <- performance(pred.B, measure='tpr', x.measure='fpr')

### SIMPLE FEATURES + RANDOM INTEGERS (CONTROL)

randnums <- sample(x=1:nrow(data))
data$randnum <- randnums

# 75/25 train/test split
index <- sample(x = 1:nrow(data), size = round(0.75 * nrow(data))) # indices of observations to be included in train set, drawing w/o replacement by default
train <- data[index,]
test <- data[-index,] # rows that train didn't include

# logistic regression
model <- glm(ie ~ gc + nstart + maxnstop + randnum, family = 'binomial', train)
p <- predict(model, test, type='response')
summary(p)
p_class <- ifelse(p > .50, 1, 0)
confusionMatrix(p_class, test[['ie']])
colAUC(p, test, plotROC = TRUE)
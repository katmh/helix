library(caret)
library(ROCR)

set.seed(500)

setwd('C:/Users/Jan_Huang/Desktop/helix_grindtime')
data <- read.csv('dataset.csv')

# add disparity and random numbers columns
maxdisps <- as.matrix(read.delim('maxdisparities.txt', header=FALSE))
data$maxdisp <- maxdisps
data <- data[!(is.na(data$maxdisp)), ] # remove rows with NA in maxdisp column
#write.csv(data, 'dataset_disp.csv'); end up with 1981 observations

randnums <- sample(x=1:nrow(data))
data$randnum <- randnums
#write.csv(data, 'dataset_disp_rand.csv')

# 75/25 train/test split
index <- sample(x = 1:nrow(data), size = round(0.75 * nrow(data))) # indices of observations to be included in train set, drawing w/o replacement by default
train <- data[index,] # 1486 observations
test <- data[-index,] # rows that train didn't include; 495 observations

# logistic regression
model <- glm(ie ~ gc + nstart + maxnstop + maxdisp + randnum, family = 'binomial', train)
p <- predict(model, test, type='response')
summary(p)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.02191 0.64747 0.48415 0.78517 0.95347
p_class <- ifelse(p > .50, 1, 0) # rn accuracy + sensitivity improve if condition is p>.55
confusionMatrix(p_class, test[['ie']])

# ROC for entire model
pred.A <- prediction(p, test$ie)
ROC.A <- performance(pred.A, measure='tpr', x.measure='fpr')
plot(ROC.A, main='ROC')

### FEATURE SELECTION

# remove redundant features
(correlationMatrix <- cor(data[,3:7]))
(highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.5))
# nstart and maxnstop are very highly correlated (0.99); generally should remove features >0.75 correlated

# rank featurs by importance, learning vector quantization (LVQ)
importance <- varImp(model, scale=FALSE)
print(importance) # maxdisp is least important omg
#            Overall
#gc        3.0367897
#nstart    4.2514813
#maxnstop 10.4119155
#maxdisp   0.1873014
#randnum   1.8962457

# recursive feature elimination (RFE), which uses random forests cross-validation
# this part is a work in progress, haven't fully interpreted it yet
control <- rfeControl(functions=rfFuncs, method='cv', number=10)
results <- rfe(data[,3:7], data[,2], sizes=c(3:7), rfeControl=control) # got 31 warnings about <=5 unique values
print(results)
predictors(results)
plot(results, type=c('g', 'o'))

### RANDOM FORESTS
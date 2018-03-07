library(caret)

setwd('C:/Users/Jan_Huang/Desktop/helix_grindtime/by_secondary_structure')
train <- read.csv('data_train_proportional.csv')
test <- read.csv('data_test_proportional.csv')
train <- subset(train, select=c('seq', 'code', 'disp', 'GC', 'ss', 'maxORF', 'ORFcover', 'maxnstop'))
data$code <- as.factor(data$code)

set.seed(500)

# np = non-proportional
index <- sample(x = 1:nrow(data), size = round(0.75 * nrow(data)))
np_train <- data[index,]
np_test <- data[-index,]

model <- glm(code ~ GC + disp, family = 'binomial', np_train)
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

a75_index <- sample(x = 1:nrow(data[data$ss=='a',]), size = round(0.75 * nrow(data[data$ss=='a',])))
a75 <- data[a75_index,]
model <- glm(code ~ GC + disp, family = 'binomial', a75)

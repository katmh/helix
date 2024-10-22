library(dplyr)
library(caret)
library(ROCR)
library(ggplot2)
library(plotly)
library(corrplot)

### QUICK START
setwd('C:/Users/Jan_Huang/Desktop/helix_grindtime/by_secondary_structure')
data <- read.csv('data_full.csv')
data <- subset(data, select=c('seq', 'code', 'disp', 'GC', 'ss', 'maxORF', 'ORFcover', 'maxnstop'))
data$code <- as.factor(data$code)

# make results reproducible
set.seed(500)

### INITIAL DATASET CLEANING
# filtering by length is done in Python
# make datasets the same size by randomly choosing sequences from bigger dataset
fraction <- nrow(coding) / nrow(noncoding)
nc_index <- sample(x = 1:nrow(noncoding), size = round(fraction * nrow(noncoding)))
noncoding <- noncoding[nc_index, ]
# remove duplicate sequences
data <- data[!duplicated(data$seq), ]

### FEATURE EXTRACTION
# values imported from files generated by Python algorithms
disp <- as.matrix(read.csv('data_full_disp.csv'))
for (i in 1:nrow(disp)) {
  data$disp[i] <- disp[i]
}
GC <- as.matrix(read.csv('data_full_GC.csv'))
for (i in 1:nrow(GC)) {
  data$GC[i] <- GC[i]
}
maxORF <- as.matrix(read.csv('data_full_maxORF.csv'))
for (i in 1:nrow(maxORF)) {
  data$maxORF[i] <- maxORF[i]
}
ORFcover<- as.matrix(read.csv('data_full_ORFcover.csv'))
for (i in 1:nrow(ORFcover)) {
  data$ORFcover[i] <- ORFcover[i]
}
maxNStop <- as.matrix(read.csv('data_full_maxNStop.csv'))
for (i in 1:nrow(maxNStop)) {
  data$maxnstop[i] <- maxNStop[i]
}
write.csv(data, 'data_full.csv')

### DESCRIPTIVE STATISTICS FOR FEATURES

# compare sequence lengths
for (i in 1:nrow(coding)) {
  coding$length[i] <- nchar(as.character(coding$seq[i]))
}
for (i in 1:nrow(noncoding)) {
  noncoding$length[i] <- nchar(as.character(noncoding$seq[i]))
}

hist(coding$length, main='Length of Coding Sequences', xlim=c(0, 100000), breaks=100, xlab='Length')
hist(noncoding$length, main='Length of Noncoding Sequences', xlim=c(0, 100000), breaks=100, xlab='Length')

# correlation matrix
feature_cols <- c(3,4,6,7,8)
M <- cor(data[feature_cols])
corrplot(M, type='upper', method='color', tl.col='black', tl.srt=45, addCoef.col='black', diag=FALSE, mar=c(0,0,1,0))

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
# conclusion: reject null that distributions are normal

fligner.test(data[data$code=='1',]$GC, data[data$code=='0',]$GC) # p=.405, variances are same
fligner.test(data[data$code=='1',]$disp, data[data$code=='0',]$disp) # p=.492, variances are same
# conclusion: variances are same

# Wilcoxon Rank Sum Test (compare mean of 2 samples)
# alternative: true location shift is not equal to 0
wilcox.test(data[data$code=='1',]$GC, data[data$code=='0',]$GC) # p=9.46e-14
wilcox.test(data[data$code=='1',]$disp, data[data$code=='0',]$disp) # p=2.2e-16
# conclusion: from diff distributions
a_b <- wilcox.test(data[data$ss=='a',]$disp, data[data$ss=='b',]$disp) # 3.460092e-76
a_ab <- wilcox.test(data[data$ss=='a',]$disp, data[data$ss=='ab',]$disp) # 3.979192e-05
a_fSS <- wilcox.test(data[data$ss=='a',]$disp, data[data$ss=='fSS',]$disp) # 0.1143285 accept null!
b_ab <- wilcox.test(data[data$ss=='b',]$disp, data[data$ss=='ab',]$disp) # 4.563859e-58
b_fSS <- wilcox.test(data[data$ss=='b',]$disp, data[data$ss=='fSS',]$disp) # 1.117752e-05
ab_fSS <- wilcox.test(data[data$ss=='ab',]$disp, data[data$ss=='fSS',]$disp) # 0.5295886 accept null!

ggplot(data=data, aes(data$disp, fill=code)) + geom_histogram(binwidth = .05) + labs(title='Histogram of Median Disparity', x='Median Disparity', y='Count') + scale_fill_discrete(name = '', labels = c('Noncoding', 'Coding'))

# by SS!
ggplot(data=data, aes(x=data$disp, y=..density.., color=code, fill=ss)) +
  xlim(0.3,1.2) +
  scale_fill_discrete(name = 'Secondary Structure', breaks = c('', 'a', 'ab', 'b', 'fSS'), labels = c('Noncoding', 'Mainly Alpha', 'Alpha Beta', 'Mainly Beta', 'Few Secondary Structures')) +
  labs(title='Median Disparity by Secondary Structure', x='Median Disparity', y='Density') +
  geom_density(alpha=0.1)

ggplot(data=data, aes(x=data$GC, y=..density.., color=code, fill=code)) +
  scale_fill_discrete(name = '', labels = c('Noncoding', 'Coding')) +
  labs(title='GC Percentage vs Coding Status', x='GC%', y='Density') +
  geom_density(alpha=0.1)

plot_ly(data=data, x=~GC, y=~disp, z=~maxORF, color=~code)

ggplot(data=data, aes(x=data$disp, y=..density.., color=ss, fill=ss)) + labs(title='Median Disparity', x='Median Disparity', y='Density') + geom_density(position='fill')

ggplot(data=data, aes(data$GC, fill=code)) + geom_histogram(binwidth = 1.15) + labs(title='Histogram of GC Percentage', x='GC Percentage', y='Count') + scale_fill_discrete(name = '', labels = c('Noncoding', 'Coding'))

ggplot(data=data, aes(data$maxORF, fill=code)) + geom_histogram() + labs(title='Histogram of Max ORF Length', x='Max ORF Length', y='Count') + scale_fill_discrete(name = '', labels = c('Noncoding', 'Coding'))

# visualize correlations btwn/among features
ggplot(data=data, aes(x=data$GC, y=data$disp)) + geom_point(aes(GC, disp, color=code)) + labs(title='Scatterplot of GC Percentage vs Median Disparity')

ggplot(data=data, aes(x=data$maxORF, y=data$disp)) + geom_point(aes(GC, disp, color=code)) + labs(title='Scatterplot')

### 75/25 TRAIN/TEST PROPORTIONAL SPLIT

# proportional split
train_indices_a <- sample(x = as.vector(which(data$ss=='a')), size = round(0.75 * nrow(data[data$ss=='a', ])))
train_indices_b <- sample(x = as.vector(which(data$ss=='b')), size = round(0.75 * nrow(data[data$ss=='b', ])))
train_indices_ab <- sample(x = as.vector(which(data$ss=='ab')), size = round(0.75 * nrow(data[data$ss=='ab', ])))
train_indices_fSS <- sample(x = as.vector(which(data$ss=='fSS')), size = round(0.75 * nrow(data[data$ss=='fSS', ])))
train_indices_nc <- sample(x = as.vector(which(data$code==0)), size = round(0.75 * nrow(data[data$code==0, ])))

train_indices_all <- c(train_indices_a, train_indices_b, train_indices_ab, train_indices_fSS, train_indices_nc)
train <- data[train_indices_all,]
test <- data[-train_indices_all,]

# save just in case :)
write.csv(train, 'data_train_proportional.csv')
write.csv(test, 'data_test_proportional.csv')

### CUTOFF OPTIMIZATION

# control: plot sensitivity and specificity vs probability threshold
# upon inspection, optimal cutoff <.5
control_pred <- prediction(control_p, test$code)
control_sens <- performance(control_pred, measure='sens', x.measure='cutoff') # true positive rate
control_spec <- performance(control_pred, measure='spec', x.measure='cutoff') # true negative rate
plot(control_sens, main='Control Model: Sensitivity and Specificity vs Cutoff', col='red', ylab='Sensitivity/Specificity')
plot(control_spec, add=TRUE, col='blue')
legend(.015, .85, legend=c("Sensitivity", "Specificity"),
       col=c("red", "blue"), lty=1:1, text.font=1)

# optimal probability threshold via largest sum
control_best.sum <- which.max(control_sens@y.values[[1]]+control_spec@y.values[[1]])
control_sens@x.values[[1]][control_best.sum] # 41.76

# optimal probability threshold via closest intersection
control_both.eq <- which.min(abs(control_sens@y.values[[1]]-control_spec@y.values[[1]]))
control_sens@x.values[[1]][control_both.eq] # 45.87

# experimental: plot sensitivity and specificity vs probability threshold
# upon inspection, optimal cutoff <.5
exp_pred <- prediction(exp_p, test$code)
exp_sens <- performance(exp_pred, measure='sens', x.measure='cutoff') # true positive rate
exp_spec <- performance(exp_pred, measure='spec', x.measure='cutoff') # true negative rate
plot(exp_sens, main='Experimental Model: Sensitivity and Specificity vs Cutoff', col='red', ylab='Sensitivity/Specificity')
plot(exp_spec, add=TRUE, col='blue')

# optimal probability threshold via largest sum
exp_best.sum <- which.max(exp_sens@y.values[[1]]+exp_spec@y.values[[1]])
exp_sens@x.values[[1]][exp_best.sum] # 39.65

# optimal probability threshold via closest intersection
exp_both.eq <- which.min(abs(exp_sens@y.values[[1]]-exp_spec@y.values[[1]]))
exp_sens@x.values[[1]][exp_both.eq] # 46.99

### LOGISTIC REGRESSION

# I'll use cutoffs optimized by largest sum

# control: standard features only (i.e. no disparity)
control_model <- glm(code ~ GC + maxORF + ORFcover + maxnstop, family = 'binomial', train)
summary(control_model) # AIC: 7644.6
control_p <- predict(control_model, test, type = 'response')
control_p_class <- ifelse(control_p > .4176, 1, 0)
confusionMatrix(control_p_class, test[['code']])
### USING .50 CUTOFF:
#          Reference
#Prediction   0   1
#         0 845 400
#         1 281 727
# Accuracy: 69.77
# Kappa: 39.55
# Sensitivity: 75.04
# Specificity: 64.51
### USING OPTIMIZED (LARGEST SUM) 41.76 CUTOFF:
#          Reference
#Prediction   0   1
#         0 727 255
#         1 399 872
# Accuracy: 70.97
# Kappa: 41.94
# Sensitivity: 64.56
# Specificity: 77.37
### USING OPTIMIZED (CLOSEST INTERSECTION) 45.87 CUTOFF:
#          Reference
#Prediction   0   1
#         0 785 341
#         1 341 786
# Accuracy: 69.73
# Kappa: 39.46
# Sensitivity: 69.72
# Specificity: 69.74

# experimental: standard features + disparity
exp_model <- glm(code ~ GC + maxORF + ORFcover + maxnstop + disp, family = 'binomial', train)
summary(exp_model) # AIC: 6493.5
exp_p <- predict(exp_model, test, type = 'response')
exp_p_class <- ifelse(exp_p > .3965, 1, 0)
confusionMatrix(exp_p_class, test[['code']])
### USING .50 CUTOFF:
#          Reference
#Prediction   0   1
#         0 920 255
#         1 206 872
# Accuracy: 79.54
# Kappa: 54.64
# Sensitivity: 81.71
# Specificity: 77.37
### USING OPTIMIZED (LARGEST SUM) 39.65 CUTOFF:
#          Reference
#Prediction   0   1
#         0 814 187
#         1 312 940
# Accuracy: 77.85
# Kappa: 55.70
# Sensitivity: 72.29
# Specificity: 83.41
### USING OPTIMIZED (CLOSEST INTERSECTION) 41.76 CUTOFF:
#          Reference
#Prediction   0   1
#         0 829 213
#         1 297 914
# Accuracy: 77.36
# Kappa: 54.73
# Sensitivity: 73.62
# Specificity: 81.10

results <- c('Control', 'Experimental')
results <- as.data.frame(results)
results$accuracy <- c(70.97, 77.85)

### TEST BY SEC STRUC

# a proportion: 0.361517976
test_a <- test[test$ss=='a', ] # 407 alpha, 407 nc
nc_indices_407 <- sample(x = as.vector(which(test$code==0)), size = 407)
test_a <- rbind(test_a, test[nc_indes_407, ])
control_p_a <- predict(control_model, test_a, type = 'response')
control_p_class_a <- ifelse(control_p_a > .4176, 1, 0)
confusionMatrix(control_p_class_a, test_a[['code']])
# Accuracy: 70.02, Sens: 66.83, Spec: 73.22
# NIR: 0.5, P-Value (Acc>NIR): <2e-16
exp_p_a <- predict(exp_model, test_a, type = 'response')
exp_p_class_a <- ifelse(exp_p_a > .3965, 1, 0)
confusionMatrix(exp_p_class_a, test_a[['code']])
# Accuracy: 78.75, Sens: 72.73, Spec: 84.77
# NIR: 0.5, P-Value (Acc>NIR): <2.2-16

# redo these bc they're messy and small mistakes were made lol
test_b <- test[test$ss=='b', ] # 323 alpha, 323 nc
nc_indices_323 <- sample(x = as.vector(which(test$code==0)), size = 323)
test_b <- rbind(test_b, test[nc_indices_323, ])
control_p_b <- predict(control_model, test_b, type = 'response')
control_p_class_b <- ifelse(control_p_b > .4176, 1, 0)
confusionMatrix(control_p_class_b, test_b[['code']])
# Accuracy: 74.3, Sens: 67.95, Spec: 73.22
# NIR: 0.642, P-Value (Acc>NIR): 3.500e-05
exp_p_b <- predict(exp_model, test_b, type = 'response')
exp_p_class_b <- ifelse(exp_p_b > .3965, 1, 0)
confusionMatrix(exp_p_class_b, test_b[['code']])
# Accuracy: 78.01, Sens: 74.25, Spec: 84.77
# NIR: .642, P-Value (Acc>NIR): <2.2e-16

test_ab <- test[test$ss=='ab', ] # 379 alpha, 379 nc
nc_indices_379 <- sample(x = as.vector(which(test$code==0)), size = 379)
test_b <- rbind(test_a, test[nc_indices_379, ])
control_p_b <- predict(control_model, test_b, type = 'response')
control_p_class_b <- ifelse(control_p_b > .4176, 1, 0)
confusionMatrix(control_p_class_b, test_b[['code']])
# Accuracy: 69.83, Sens: 67.95, Spec: 73.22
# NIR: 0.642, P-Value (Acc>NIR): 3.500e-05
exp_p_b <- predict(exp_model, test_b, type = 'response')
exp_p_class_b <- ifelse(exp_p_b > .3965, 1, 0)
confusionMatrix(exp_p_class_b, test_b[['code']])
# Accuracy: 78.01, Sens: 74.25, Spec: 84.77
# NIR: .642, P-Value (Acc>NIR): <2.2e-16

### ROC CURVES
control_roc <- performance(control_pred, measure = 'tpr', x.measure = 'fpr')
plot(control_roc, main='Control: ROC Curve')
abline(a=0,b=1)
control_auc <- performance(control_pred, measure = 'auc')
control_auc@y.values # 0.7895

exp_roc <- performance(exp_pred, measure = 'tpr', x.measure = 'fpr')
plot(exp_roc, main='ROC Curve', col = 'red')
plot(add = TRUE, control_roc, col = 'blue')
abline(a=0,b=1)
legend(0.35, .2, legend=c("Control (AUROC: 0.7895)", "Experimental (AUROC: 0.8605)"),
       col=c("blue", "red"), lty=1:1, cex=0.8)
exp_auc <- performance(exp_pred, measure = 'auc')
exp_auc@y.values # 0.8605

### RANDOM NOISE
# if any features are less important than the random noise variable, they are insignificant

# generate random numbers
randnums <- sample(x=1:nrow(data), size=nrow(data))
for (i in 1:nrow(train)) {
  train$randnum[i] <- randnums[i]
}
for (i in 1:nrow(test)) {
  test$randnum[i] <- randnums[i+nrow(train)]
}

# train model
noise_model <- glm(code ~ GC + maxORF + ORFcover + maxnstop + disp + randnum, family = 'binomial', data=train)
noise_p <- predict(noise_model, test, type = 'response')

# optimize cutoff
noise_pred <- prediction(noise_p, test$code)
noise_sens <- performance(noise_pred, measure='sens', x.measure='cutoff') # true positive rate
noise_spec <- performance(noise_pred, measure='spec', x.measure='cutoff') # true negative rate
noise_best.sum <- which.max(noise_sens@y.values[[1]]+noise_spec@y.values[[1]])
noise_sens@x.values[[1]][exp_best.sum] # 39.93

# test model
noise_p_class <- ifelse(noise_p > .3993, 1, 0)
confusionMatrix(noise_p_class, test[['code']])
#          Reference
#Prediction   0   1
#         0 816 187
#         1 310 940
# Accuracy: 77.94
# Kappa: 55.88
# Sensitivity: 72.47
# Specificity: 83.41
(noise_varImp <- varImp(noise_model, scale=FALSE))
#GC       14.3785158
#maxORF   18.3890397
#ORFcover 14.7707025
#maxnstop 11.0902717
#disp     27.9960533
#randnum   0.3756356

### RANK FEATURES BY IMPORTANCE
# if feature less important than random noise, then it's not significant
# uses learning vector quantization (LVQ)
(control_varImp <- varImp(control_model, scale=FALSE))
#GC       10.77046
#maxORF   18.14501
#ORFcover 17.09556
#maxnstop 12.15270
control_varImp <- as.data.frame(control_varImp$Overall)
control_varImp$Feature <- c('GC Percentage', 'Maximum ORF Length', 'ORF Coverage', 'Maximum Occurrences of Stop Codons')
colnames(control_varImp) <- c('Importance', 'Feature')
ggplot(control_varImp, aes(x = reorder(Feature, -Importance), y = Importance, fill = Feature)) +
  geom_col() +
  scale_x_discrete(labels = function(x) lapply(strwrap(x, width = 17.5, simplify = FALSE), paste, collapse="\n")) +
  scale_y_continuous(limits=c(0,30)) +
  labs(x = 'Feature') +
  geom_text(aes(label=substr(as.character(Importance),1,5)), position=position_dodge(width=0.9), vjust=-1)
(exp_varImp <- varImp(exp_model, scale=FALSE))
#GC       14.37687
#maxORF   18.38776
#ORFcover 14.77688
#maxnstop 11.09030
#disp     27.99592
exp_varImp <- as.data.frame(exp_varImp$Overall)
exp_varImp$Feature <- c('GC Percentage', 'Maximum ORF Length', 'ORF Coverage', 'Maximum Occurrences of Stop Codons', 'Median Disparity')
colnames(exp_varImp) <- c('Importance', 'Feature')
ggplot(exp_varImp, aes(x = reorder(Feature, -Importance), y = Importance, fill = Feature)) +
  geom_col() +
  scale_x_discrete(labels = function(x) lapply(strwrap(x, width = 15, simplify = FALSE), paste, collapse="\n")) +
  scale_y_continuous(limits=c(0,30)) +
  labs(x = 'Feature') +
  geom_text(aes(label=substr(as.character(Importance),1,5)), position=position_dodge(width=0.9), vjust=-1)

### CROSS VALIDATION
# 10-fold cV on control, using Bayes GLM
(control_cv <- train(
  code ~ GC + maxORF + ORFcover + maxnstop, data,
  method = 'glm',
  trControl = trainControl(
    method = 'cv', number = 10,
    verboseIter = TRUE
  )
))
# accuracy: 70.46
# kappa: 40.92 - borderline poor
# Kappa = (observed accuracy - expected accuracy)/(1 - expected accuracy)

# 10-fold cV on control, using Bayes GLM
(exp_cv <- train(
  code ~ GC + maxORF + ORFcover + maxnstop + disp, data,
  method = 'glm',
  trControl = trainControl(
    method = 'cv', number = 10,
    verboseIter = TRUE
  )
))
# accuracy: 78.96
# kappa: 57.92 - moderate/good

### ANOVA
# all features add significant improvement to the model
two_anova <- anova(control_model, exp_model, test = 'Chisq')
two_anova$`Pr(>Chi)` # 3.894958e-254

control_anova <- anova(control_model, test = 'Chisq')
control_anova$`Pr(>Chi)` # GC: 6.671715e-91, maxORF: 4.640148e-155, ORFcover: 1.872042e-85, maxnstop: 3.062351e-48

exp_anova <- anova(exp_model, test = 'Chisq')
exp_anova$`Pr(>Chi)` # GC: 6.671715e-91, maxORF: 4.640148e-155, ORFcover: 1.872042e-85, maxnstop: 3.062351e-48, disp: 3.894958e-254

mean(data[data$ss=='a',]$disp) # .730
mean(data[data$ss=='b',]$disp) # .658
mean(data[data$ss=='ab',]$disp) # .707
mean(data[data$ss=='fSS',]$disp) # .712
mean(data[data$code==0,]$disp) # .615
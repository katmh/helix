data <- read.table('statistics.txt', sep = '\n', skip = 4)
data <- as.numeric(unlist(data))

hist(data,
     breaks = 100,
     xlab = "Maximum Values",
     main = "Histogram of Maximum Disparity Values",
     xlim = c(0.5,2.5))

mean <- mean(data)
mean <- substr(as.character(mean), 0, 5)

stdev <- sd(data)
stdev <- substr(as.character(stdev), 0, 5)

text(x=200, labels=paste('Mean: ', mean, '\nSD: ', stdev))

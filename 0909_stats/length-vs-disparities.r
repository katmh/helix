data <- read.csv('length-vs-disparities.txt', header = TRUE)

plot(data$length, data$maximum)

plot(x = data$length,
     y = data$maximum,
     main = 'Sequence Length vs. Maximum Disparity Value',
     xlab = 'Sequence Length (Number of Residues)',
     ylab = 'Maximum Disparity Value',
     ylim = c(0, 2.5))
text(55, .1, labels=paste('Mean: ', mean(data$maximum), '\nSD: ', sd(data$maximum)))

plot(x = data$length,
     y = data$mean,
     main = 'Sequence Length vs. Mean Disparity Value',
     xlab = 'Sequence Length (Number of Residues)',
     ylab = 'Mean Disparity Value',
     ylim = c(0, 2.5))
text(55, .1, labels=paste('Mean: ', mean(data$mean), '\nSD: ', sd(data$mean)))

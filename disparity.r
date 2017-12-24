setwd('~/helix2/helix')

mm_mainly_alpha <- read.csv('mm_mainly_alpha.csv', header=TRUE, sep=',')
mm_mainly_beta <- read.csv('mm_mainly_beta.csv', header=TRUE, sep=',')
mm_alpha_beta <- read.csv('mm_alpha_beta.csv', header=TRUE, sep=',')
mm_few_ss <- read.csv('mm_few_ss.csv', header=TRUE, sep=',')
mm_NONCODE <- read.csv('mm_NONCODE.csv', header=TRUE, sep=',')

# Mean and Maximum for Mainly Alpha
hist(mm_mainly_alpha$mean,
     breaks = 100,
     xlim = c(0,2.75),
     col = 'blue',
     main = paste('Mean and Maximum Disparities in Mainly-Alpha Protein Chains'),
     xlab = 'Disparity Value')
hist(mm_mainly_alpha$max,
     breaks = 100,
     col = 'red',
     add = TRUE)
legend('topright',
       inset = .1,
       c('Mean','Maximum'),
       fill = c('blue', 'red'))

# Mean and Maximum for Mainly Beta
hist(mm_mainly_beta$mean,
     breaks = 100,
     xlim = c(0,2.75),
     col = 'blue',
     main = paste('Mean and Maximum Disparities in Mainly-Beta Protein Chains'),
     xlab = 'Disparity Value')
hist(mm_mainly_beta$max,
     breaks = 100,
     col = 'red',
     add = TRUE)
legend('topright',
       inset = .1,
       c('Mean','Maximum'),
       fill = c('blue', 'red'))

# Mean and Maximum for Mainly Beta
hist(mm_alpha_beta$mean,
     breaks = 100,
     xlim = c(0,2.75),
     col = 'blue',
     main = paste('Mean and Maximum Disparities in Alpha-Beta Protein Chains'),
     xlab = 'Disparity Value')
hist(mm_alpha_beta$max,
     breaks = 100,
     col = 'red',
     add = TRUE)
legend('topright',
       inset = .1,
       c('Mean','Maximum'),
       fill = c('blue', 'red'))

# Mean and Maximum for Few Secondary Structures
hist(mm_few_ss$mean,
     breaks = 100,
     xlim = c(0,2.75),
     col = 'blue',
     main = paste('Mean and Maximum Disparities in Protein Chains with Few Secondary Structures'),
     xlab = 'Disparity Value')
hist(mm_few_ss$max,
     breaks = 100,
     col = 'red',
     add = TRUE)
legend('topright',
       inset = .1,
       c('Mean','Maximum'),
       fill = c('blue', 'red'))

# Mean and Maximum for Non-Coding RNA
hist(mm_NONCODE$mean,
     breaks = 100,
     xlim = c(0,2.75),
     col = 'blue',
     main = paste('Mean and Maximum Disparities: Non-Coding RNA'),
     xlab = 'Disparity Value')
hist(mm_NONCODE$max,
     breaks = 100,
     col = 'red',
     add = TRUE)
legend('topright',
       inset = .1,
       c('Mean','Maximum'),
       fill = c('blue', 'red'))

# Maximum for Mainly Alpha, Mainly Beta, and Alpha-Beta
hist(mm_alpha_beta$max,
     breaks = 100,
     xlim = c(0,3),
     col = rgb(0,1,0,.5),
     main = paste('Maximum Disparities of Mainly-Alpha, Mainly-Beta, and Alpha-Beta Protein Chains'),
     xlab = 'Disparity Value')
hist(mm_mainly_beta$max,
     breaks = 100,
     col = rgb(0,0,1,.5),
     add = TRUE)
hist(mm_mainly_alpha$max,
     breaks = 100,
     col = rgb(1,0,0,.5),
     add = TRUE)
legend('topright',
       inset = .1,
       c('Mainly Beta','Mainly Alpha', 'Alpha Beta'),
       fill = c(rgb(0,0,1,.5), rgb(1,0,0,.5), rgb(0,1,0,.5)))

# Mean for Mainly Alpha, Mainly Beta, and Alpha-Beta
hist(mm_alpha_beta$mean,
     breaks = 100,
     xlim = c(0,2),
     col = rgb(0,1,0,.75),
     main = paste('Mean Disparities: Mainly-Alpha, Mainly-Beta, and Alpha-Beta Protein Chains'),
     xlab = 'Disparity Value')
hist(mm_mainly_alpha$mean,
     breaks = 100,
     col = rgb(1,0,0,.75),
     add = TRUE)
hist(mm_mainly_beta$mean,
     breaks = 100,
     col = rgb(0,0,1,.75),
     add = TRUE)
legend('topright',
       inset = .1,
       c('Mainly Beta','Mainly Alpha', 'Alpha Beta'),
       fill = c(rgb(0,0,1,.75), rgb(1,0,0,.75), rgb(0,1,0,.75)))

# Maximum for Alpha-Beta and Few Secondary Structures
hist(mm_alpha_beta$max,
     breaks = 100,
     xlim = c(0,3),
     col = rgb(0,1,0,.5),
     main = paste('Maximum Disparities of Alpha-Beta and Few-Secondary-Structure Protein Chains'),
     xlab = 'Disparity Value')
hist(mm_few_ss$max,
     breaks = 100,
     col = rgb(1,0,0,.5),
     add = TRUE)
legend('topright',
       inset = .05,
       c('Alpha Beta', 'Few Secondary Structures'),
       fill = c(rgb(0,1,0,.5), rgb(1,0,0,.5)))

# Maximum for Alpha-Beta, Few Secondary Structures, and Non-Coding RNA
hist(mm_alpha_beta$max,
     breaks = 100,
     xlim = c(0,3),
     col = rgb(0,1,0,.75),
     main = paste('Maximum Disparities: Alpha-Beta, Few Secondary Structures, and Non-Coding RNA'),
     xlab = 'Disparity Value')
hist(mm_NONCODE$max,
     breaks = 100,
     col = rgb(0,0,1,.75),
     add = TRUE)
hist(mm_few_ss$max,
     breaks = 100,
     col = rgb(1,0,0,.75),
     add = TRUE)
legend('topright',
       inset = .05,
       c('Alpha Beta', 'Few Secondary Structures', 'Non-Coding RNA'),
       fill = c(rgb(0,1,0,.75), rgb(1,0,0,.75), rgb(0,0,1,.75)))

# Means for Alpha-Beta, Few Secondary Structures, and Non-Coding RNA
hist(mm_alpha_beta$mean,
     breaks = 100,
     xlim = c(0,1.5),
     col = rgb(0,1,0,.75),
     main = paste('Mean Disparities: Alpha-Beta, Few Secondary Structures, and Non-Coding RNA'),
     xlab = 'Disparity Value')
hist(mm_NONCODE$mean,
     breaks = 100,
     col = rgb(0,0,1,.75),
     add = TRUE)
hist(mm_few_ss$mean,
     breaks = 100,
     col = rgb(1,0,0,.75),
     add = TRUE)
legend('topright',
       inset = .05,
       c('Alpha Beta', 'Few Secondary Structures', 'Non-Coding RNA'),
       fill = c(rgb(0,1,0,.75), rgb(1,0,0,.75), rgb(0,0,1,.75)))
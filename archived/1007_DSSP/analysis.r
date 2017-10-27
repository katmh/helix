ID1ajqB01 <- read.csv('1ajqB01.csv', header=TRUE, sep=',')

plot(ID1ajqB01$ss, ID1ajqB01$disp)

plot(x = ID1ajqB01$ss,
     y = ID1ajqB01$disp,
     xlab = 'Secondary Structure',
     ylab = 'Disparity Value',
     main = 'Disparity Levels for Different Secondary Structures')

# isolate residues categorized as helices
ID1ajqB01.H <- subset(ID1ajqB01, ID1ajqB01$ss == 'H')
# isolate residues categorized as beta bridge
ID1ajqB01.B <- subset(ID1ajqB01, ID1ajqB01$ss == 'B')

hist(ID1ajqB01.H$disp)

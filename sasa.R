setwd('~/helix2/helix')

pp_beta <- read.csv('pp_mainly_beta.csv')
pp_alpha <- read.csv('pp_mainly_alpha.csv')
pp_ab <- read.csv('pp_ab.csv')
sasa_beta <- read.csv('sasa_mainly_beta.csv')
sasa_alpha <- read.csv('sasa_mainly_alpha.csv')
sasa_ab <- read.csv('sasa_ab.csv')

# merge data frames
pps_beta <- cbind(pp_beta, sasa_beta, by = 'ID')
pps_alpha <- cbind(pp_alpha, sasa_alpha, by = 'ID')
# remove redundant ID column?

scatter.smooth(x = pps_beta$percentPhilic,
               y = pps_beta$normSASA,
               xlab = 'Percent Hydrophilic Residues',
               ylab = 'Normalized Solvent-Accessible Surface Area',
               main = '% Hydrophilicity vs. Accessible Surface Area: Mainly Beta Proteins')

scatter.smooth(x = pps_beta$percentPhobic,
               y = pps_beta$normSASA,
               xlab = 'Percent Hydrophobic Residues',
               ylab = 'Normalized Solvent-Accessible Surface Area',
               main = 'Hydrophobicity Versus Accessible Surface Area: Mainly Beta Proteins')

scatter.smooth(x = pps_alpha$percentPhobic,
               y = pps_alpha$normSASA,
               xlab = 'Percent Hydrophobic Residues',
               ylab = 'Normalized Solvent-Accessible Surface Area',
               main = 'Hydrophobicity Versus Accessible Surface Area: Mainly Alpha Proteins')

hist(x = pp_beta$percentPhobic,
     xlab = 'Percent Hydrophobicity',
     breaks = 100,
     main = 'Histogram of Percent Hydrophobicity in Mainly Beta Proteins')

hist(x = pps_beta$normSASA,
     xlab = 'Normalized Solvent-Accessible Surface Area',
     breaks = 100,
     xlim = c(2,12),
     main = 'Histogram of Normalized SASA in Mainly Beta Proteins')

mean(pp_beta$percentPhobic) == 1 - mean(pp_beta$percentPhilic) # true

sd(pp_beta$percentPhobic)
sd(pp_beta$percentPhilic)

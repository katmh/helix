# Calculate hexamer frequency from coding and noncoding training data sets

hexFreq = {}

hexamers = []

nucleotides = ['A', 'T', 'C', 'G']



def myFunc(i):
	with open(i) as f:
		for line in f:

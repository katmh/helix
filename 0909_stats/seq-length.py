# originally looked at outliers (2fmlA01 max: 0.86, 4ad9A01 max: 2.25)
# 2fmlA01 is much shorter, so I want to analyze seq length vs max/mean disparity value

domains = []
lengths = []
maxima = []
means = []

def mean(numbers):
	return float( sum(numbers) / len(numbers) )

with open('../0904_dataset_threshold/disparities-n68.txt') as f:
	
	for i, line in enumerate(f):

		if i % 2 == 0: # domains
			domains.append(line[:7])
			maxima.append(line[13:-1])
	
		if i % 2 == 1:
			values = []

			for i in line[:-2].split(' '):
				values.append(float(i))

			lengths.append( len(values) + 15 ) # amino acids in sequence
			means.append( mean(values) )

# write csv
with open('length-vs-disparities.txt', 'a') as f:
	f.write('domain,length,maximum,mean\n')
	for i in range(0, len(domains) - 1):
		f.write(domains[i] + ',' + str(lengths[i]) + ',' + str(maxima[i]) + ',' + str(means[i]))
		f.write('\n')
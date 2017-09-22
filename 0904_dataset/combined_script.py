import math
import matplotlib.pyplot as plt
import statistics

def createDataset( divisor ):
	allHeader = []
	allSeq = []

	with open('cath-dataset-nonredundant-S20.fa') as f:
		for i, line in enumerate(f):
			if i % 2 == 0:
				allHeader.append(line)
			if i % 2 == 1:
				allSeq.append(line)

	smallerListHeader = []
	smallerListSeq = []

	if divisor != 1:
		for i, head in enumerate(allHeader):
			if i % divisor == 0:
				smallerListHeader.append(head)
				smallerListSeq.append(allSeq[i])
	else:
		smallerListHeader = allHeader
		smallerListSeq = allSeq

	finalListHeader = []
	finalListSeq = []
	excluded = []

	with open('cath-domain-list-S100.txt') as f:
		for line in f:
			for i, head in enumerate(smallerListHeader):
				if head[12:19] == line[0:7]:
					if (line[12:13] == '1') or (line[12:13] == '3'):
						finalListHeader.append(head)
						finalListSeq.append(smallerListSeq[i])
					else:
						excluded.append(line[0:7] + ' ' + str(line[12:13]))

	# output
	with open('exclude24.txt', 'a') as ex:
		ex.write('n = ' + str(len(excluded)) + '\n')
		for record in excluded:
			ex.write(record + '\n')

	with open('final-dataset.txt', 'a') as final:
		final.write('n = ' + str(len(finalListSeq) - 1) + '\n')
		for i in range(0, len(finalListSeq) - 1):
			final.write(finalListHeader[i][12:19] + ' ' + finalListSeq[i])

def slide( record, sectionLength ):
	# hydropathy values for each amino acid
	table = {'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9, 'A': 1.8, 'G': -0.4, 'T': -0.7, 'W': -0.9, 'S': -0.8, 'Y': -1.3, 'P': -1.6, 'H': -3.2, 'E': -3.5, 'D': -3.5, 'N': -3.5, 'Q': -3.5, 'K': -3.9, 'R': -4.5, 'X': 0}

	pdb = record[0:7]
	seq = record[8:-1]

	disparities = []

	# iterate through protein sequence to get starting resdiues of each window
	for i, slideStart in enumerate( seq[0 : len(seq) - sectionLength - 1] ):

		# variables to be populated
		hydro = []
		sumX = 0
		sumY = 0

		# from the starting residue, iterate through next n-1 residues, n=sectionLength
		for n, residue in enumerate( seq[ i : i + sectionLength - 1 ] ):

			# map residue to hydropathy value
			hydro.append(table[residue])

			# represent top-down view of helix as circle
			positionX = math.cos( math.radians(-100*n) ) # n starts at 0
			positionY = math.sin( math.radians(-100*n) )

			sumX += hydro[n] * positionX
			sumY += hydro[n] * positionY

		vectorLength = math.sqrt(sumX**2 + sumY**2) # distance formula

		disparities.append(vectorLength / float( sectionLength ))

	with open('disparities.txt', 'a') as f:
		f.write(pdb + ' max: ' + str(max(disparities)) + '\n')
		for value in disparities:
			f.write(str(value) + ' ')
		f.write('\n')

def analyze():
	maxima = []

	with open('disparities.txt') as f:
		for i, line in enumerate(f):
			if i % 2 == 0:
				maxima.append(float(line[13:-1]))

	with open('statistics.txt', 'a') as f:
		f.write('n = ' + str(len(maxima)) + '\n')
		f.write('Mean: ' + str(statistics.mean(maxima)) + '\n')
		f.write('Standard Deviation: ' + str(statistics.stdev(maxima)) + '\n\n')
		for value in maxima:
			f.write(str(value) + '\n')

createDataset(1)

proteins = []

with open('final-dataset.txt') as data:
	for i, line in enumerate(data):
		if i > 0:
			proteins.append(line)

for i, protein in enumerate(proteins):
	slide(protein, 15)

analyze()
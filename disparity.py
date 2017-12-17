import math
import statistics

# parameters: polypeptide sequence in letters, size of sliding window
def slide(seq, size):
	# Kyte and Doolittle 1982 hydropathy index
	# assign 0 to X
	# B = N or D = -3.5
	table = {'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9, 'A': 1.8, 'G': -0.4, 'T': -0.7, 'W': -0.9, 'S': -0.8, 'Y': -1.3, 'P': -1.6, 'H': -3.2, 'E': -3.5, 'D': -3.5, 'N': -3.5, 'Q': -3.5, 'K': -3.9, 'R': -4.5, 'X': 0, 'U': 2.5, 'O': -3.9, 'B': -3.5}

	disparities = []

	# iterate thru seq and get starting residues of each window
	for i, start in enumerate( seq[0 : len(seq) - (size - 1)] ):

		# variables to be populated
		hydro = []
		sumX = 0
		sumY = 0

		# from start residue, iterate thru next size-1 residues
		for n, residue in enumerate( seq[i : i + size] ):

			# map residue to hydropathy index
			hydro.append(table[residue])

			# helical wheel (top-down view of helix); residues ideally every 100deg
			# right-handed helix means -100deg, but pos/neg doesn't matter for vector magnitude
			positionX = math.cos( math.radians(100 * n) )
			positionY = math.sin( math.radians(100 * n) )

			sumX += hydro[n] * positionX
			sumY += hydro[n] * positionY

			# disparity value D is magnitude of X and Y vectors (sumX and sumY, respectively)
			disp = math.sqrt(sumX**2 + sumY**2)

		# take the avg disparity of the window and assign it to the starting residue
		disparities.append(disp / size)

	return disparities

# parameters: input CSV file w/ IDs, sequences, etc.; output file
def generateDisparities(i, o):
	with open(i) as csv:
		with open(o, 'w') as out:
			for i, line in enumerate(csv):
				if i != 0: # skip header
					if int(line.split(',')[2].strip('"')) > 20: # chain length must be greater than 20
						print(i)
						# PDB ID + chain ID
						ID = line.split(',')[0].strip('"') + line.split(',')[1].strip('"')
						out.write(ID + ',')

						seq = line.split(',')[3].split('#')[0].lstrip('"')
						for value in slide(seq, 15):
							out.write(str(value) + ',')

						out.write('\n')

# parameters: input disparity, output mean max (mm) file
def generateMeanMax(i, o):
	with open(i) as disp:
		with open(o, 'w') as out:
			out.write('ID,mean,max\n')
			for i, line in enumerate(disp):
				print(i)
				avg = statistics.mean( list(map(float, line.split(',')[1:-1])) ) # skip ID and \n char
				mx = max( line.split(',')[1:-1] )
				out.write( line.split(',')[0] + ',' + str(avg) + ',' + str(mx) + '\n' )

generateDisparities('seq_few_ss.csv', 'disp_few_ss.txt')
generateMeanMax('disp_few_ss.txt', 'mm_few_ss.csv')
import matplotlib.pyplot as plt
import numpy as np
import statistics
import math

# melittin, used as ampipathic example in Eisenberg 1982
seq = 'GIGAVLKVLTTGLPALISWIKRKRQQ'
# 26 residues

table = {'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9, 'A': 1.8, 'G': -0.4, 'T': -0.7, 'W': -0.9, 'S': -0.8, 'Y': -1.3, 'P': -1.6, 'H': -3.2, 'E': -3.5, 'D': -3.5, 'N': -3.5, 'Q': -3.5, 'K': -3.9, 'R': -4.5, 'X': 0, 'U': 2.5, 'O': -3.9, 'B': -3.5, '-': 0}

hydro = []

for residue in seq:
	hydro.append(table[residue])
'''
plt.plot(hydro)
plt.ylabel('hydropathy')
plt.xlabel('residue index')
plt.title('melittin (2mlt)')
plt.xticks(np.arange(0, len(seq), 3.6))
plt.show()

print('Mean Hydropathy: ' + str(statistics.mean(hydro)))
print(hydro)

# Eisenberg mean helical hydrophobic moment (mhhm)

# helical hydrophobic moment
hhmTerms = []

directions = []

for i, residue in enumerate(seq):
	x = math.cos( math.radians( -i * 100 ) )
	y = math.sin( math.radians( -i * 100 ) )

	# vector magnitude
	directions.append( math.sqrt( x**2 + y**2 ) )

for i, residue in enumerate(seq):
	hhmTerms.append(directions[i] * hydro[i])

print(directions)
print(hhmTerms)

plt.plot(hhmTerms)
plt.ylabel('helical hydrophobic moment')
plt.xlabel('residue index')
plt.title('melittin (2mlt) - hhm')
plt.xticks(np.arange(0, len(seq), 3.6))
plt.show()
'''
# parameters: polypeptide sequence in letters, size of sliding window
def slide(seq, size):
	# Kyte and Doolittle 1982 hydropathy index
	# assign 0 to X and stop codon
	# B = N/D = -3.5
	table = {'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9, 'A': 1.8, 'G': -0.4, 'T': -0.7, 'W': -0.9, 'S': -0.8, 'Y': -1.3, 'P': -1.6, 'H': -3.2, 'E': -3.5, 'D': -3.5, 'N': -3.5, 'Q': -3.5, 'K': -3.9, 'R': -4.5, 'X': 0, 'U': 2.5, 'O': -3.9, 'B': -3.5, '-': 0}

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

print(slide(seq, 3))
print(slide('LOOLOOLOOLOOL', 3))
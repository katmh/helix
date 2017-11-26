import math
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
rcParams['axes.titlepad'] = 10
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# parameters: input file name, divisor (get every n proteins), output file name (add, not overwrite)
def createDataset(i, n, o):
	with open(o, 'a') as output:
		with open(i) as input:
			for i, line in enumerate(input):
				# remember line 0 is header
				if i != 0 and i % n == 0:
					output.write(line)

#createDataset('alphaBeta.csv', 500, 'combinedAlphas.csv')
#createDataset('mainlyAlpha.csv', 500, 'combinedAlphas.csv')

# parameters: polypeptide sequence in letters, size of sliding window
def slide(seq, size):
	# Kyte and Doolittle 1982 hydropathy index
	# assign 0 to X
	table = {'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9, 'A': 1.8, 'G': -0.4, 'T': -0.7, 'W': -0.9, 'S': -0.8, 'Y': -1.3, 'P': -1.6, 'H': -3.2, 'E': -3.5, 'D': -3.5, 'N': -3.5, 'Q': -3.5, 'K': -3.9, 'R': -4.5, 'X': 0}

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

# parameters: input file w/ IDs and sequences, output file
def calculateDisparities(i, o):
	with open(o, 'w') as out:
		with open(i) as seqs:
			for i, line in enumerate(seqs):
				if i != 0:
					# line.split(',')[0] is PDB ID w/ surrounded by double quotes
					out.write((line.split(',')[0]).strip('"') + ',')
					
					seq = line.split(',')[2].split('#')[0].lstrip('"')
					for value in slide( seq, 15 ):
						out.write(str(value) + ',')
					
					out.write('\n')

#calculateDisparities('combinedAlphas.csv', 'disparities.txt')

# parameters: input file (CSV of PDB ID, chain ID, seqSS), output file (ID, SS)
def getSS(i, o):
	with open(o, 'w') as IDSS:
		with open(i) as CSV:
			for i, line in enumerate(CSV):
				if i != 0: # skip header
					IDSS.write(line.split(',')[0].strip('"') + ',') # write PDB ID and comma
					IDSS.write('"' + line.split('#')[1]) # write SS surrounded by double quotes, then newline

#getSS('combinedAlphas.csv', 'IDSS.csv')

# parameters: input SS; returns list
def generateSSList(SS):
	SSRegions = []

	for i, char in enumerate(SS):
		if i == 0:
			SSRegions.append((i, char))
		elif SS[i] != SS[i-1]:
			SSRegions.append((i, char))
		elif SS[i] == SS[i-1]:
			pass

	return SSRegions

def colorBySS(ID):
	# generate dictionary from SS
	with open('IDSS.csv') as IDSS:
		for line in IDSS:
			if line[0:4] == ID: # find line w/ same ID
				SSList = generateSSList(line.split(',')[1]) # arg is SS info

	for i, pair in enumerate(SSList):
		startRes = pair[0]
		if i < len(SSList) - 1:
			endRes = list(SSList)[i+1][0]
		elif i == len(SSList) - 1: # if last, region extends to the end
			endRes = len(SSList) - 1

		if pair[1] == 'H': # alpha helix
			# draw vertical span (rect) from xmin to xmax; color-coded by SS
			plt.axvspan(startRes, endRes, facecolor='m', alpha=0.4) # magenta
		if pair[1] == 'E': # beta strand/ladder
			plt.axvspan(startRes, endRes, facecolor='y', alpha=0.4) # yellow

# parameters: input file (IDs and lists of disparities), output file (PDF of plots, 1 plot per ID)
def plotDisparities(i, o):
	pdf = PdfPages(o)
	
	with open(i) as disp:

		plt.style.use('seaborn')

		for n, line in enumerate(disp):
			ID = line.split(',')[0]
			y = line.split(',')[1:-1] # last item is newline char

			fig = plt.figure(1, figsize=(10,5))

			colorBySS(ID)

			plt.subplot(111)
			plt.plot(y, c='k', linewidth=1.0)
			plt.subplots_adjust(left=0.1, right=0.95, top=0.89, bottom=0.14) # coordinate system of the figure; (0,0) is bottom left, (1,1) top right

			plt.title(ID, fontsize=14, fontweight='bold')
			plt.ylabel('Average Disparity in Window', labelpad=7.5)
			plt.xlabel('First Residue in Window', labelpad=7.5)

			pdf.savefig(fig)

			plt.clf() # clear figure to avoid overwriting

	pdf.close()

#plotDisparities('disparities.txt', 'disparities-AB-shade.pdf')
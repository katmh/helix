import math
import os
import subprocess

# Create dataset of domains containing alpha helices
# From cath-dataset-nonredundant-S20.fa, exclude classes B and D (mainly beta, few secondary structures)

# original: File path of original dataset
# domainList: File path of CATH domain list to check classifications
# output: File path of output dataset
# NOTE: this function currently takes 3 min to process CATH S20
def createDataset(original, domainList, output):

	# separate original dataset into lists of headers and sequences
	allHeaders = []
	allSeq = []

	with open(original) as originalSet:
		for i, line in enumerate(originalSet):
			if i % 2 == 0: # headers
				allHeaders.append(line[14:21])
			if i % 2 == 1:
				allSeq.append(line[:-1]) # remove newline char

	# check classification of domains (domain IDs are in headers) using domainList, then write output file

	with open(output, 'w') as outputFile:
		with open(domainList) as domains:
			for line in domains:
				for i, ID in enumerate(allHeaders):
					if line[0:7] == ID: # match IDs
						if (line[12:13] == '1') or (line[12:13] == '3'):
							# output file format: ID,seq,class,architecture
							outputFile.write(ID + ',' + allSeq[i] + ',' + str(line[12]) + ',' + str(line[17:19]) + '\n')

#createDataset('cath-dataset-nonredundant-S40-v4_2_0.fa', 'cath-domain-list-S100-v4_2_0.txt', 'dataset.txt')

# seq: polypeptide sequence in letters
# size: size of sliding window (e.g. 15 residues)
def slide( seq, size ):
	# Kyte and Doolittle 1982 hydropathy index
	table = {'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9, 'A': 1.8, 'G': -0.4, 'T': -0.7, 'W': -0.9, 'S': -0.8, 'Y': -1.3, 'P': -1.6, 'H': -3.2, 'E': -3.5, 'D': -3.5, 'N': -3.5, 'Q': -3.5, 'K': -3.9, 'R': -4.5, 'X': 0}

	disparities = []

	# iterate through sequence and get starting residues of each window
	for i, start in enumerate( seq[0 : len(seq) - (size - 1)] ):

		# variables to be populated
		hydro = []
		sumX = 0
		sumY = 0

		# from the starting residue, iterate through next size-1 residues (e.g. next 14 residues, since starting residue = 1)
		for n, residue in enumerate( seq[i : i + size] ):

			# map residue to hydropathy index
			hydro.append(table[residue])

			# helical wheel; represent top-down view of helix as circle, w/ residues occurring every 100deg
			# right-handed helix means -100deg, but we're using vectors, so doesn't matter pos/neg
			positionX = math.cos( math.radians(100 * n) )
			positionY = math.sin( math.radians(100 * n) )

			sumX += hydro[n] * positionX
			sumY += hydro[n] * positionY

			# disparity value D is magnitude (vector length) of X and Y vectors (sumX and sumY, respectively)
			disp = math.sqrt(sumX**2 + sumY**2)

		disparities.append(disp / size) # average it

	return disparities

#input: file path of dataset generated from createDataset()
# output: file path of output
# Note: function is written to input/output files w/ multiple polypeptides, not one by one (see slide function)
def calculateDisparities(input, output):
	with open(input) as dataset:
		with open(output, 'w') as out:
			for line in dataset:
				out.write(line.split(',')[0] + ',')
				for value in slide( line.split(',')[1], 15 ):
					out.write(str(value) + ',') # note: unnecessary comma at end
				out.write('\n')

#calculateDisparities('dataset.txt', 'disparities.txt')

# IDlist: list of IDs (not file)
# Note: need dompdb folder, which is very large (download from CATH)
def generateDSSP(IDlist):
	for ID in IDlist:
		subprocess.call('mkdssp -i dompdb/' + ID + ' -o dssp/' + ID + '.dssp', shell=True)

IDs = []

with open('dataset.txt') as dataset:
	for line in dataset:
		IDs.append(line.split(',')[0])

#generateDSSP(IDs)

# CSV format: header (ID); residueNum,ss,disparity
# disparity will be empty if it's last 14 residues of sequence
# ss will be ' ' if DSSP couldn't detect or it has no ss classification
# output: file path
def generateCSV(ID):
	ss = getSS('dssp/' + ID + '.dssp')
	disp = getDisp(ID, 'disparities.txt')

	with open('csv/' + ID + '.csv', 'w') as outputFile:
		outputFile.write(ID + '\n')

		for i, value in enumerate(ss):
			outputFile.write(str(i) + ',' + value)
			if i in range(0, len(disp)):
				outputFile.write(',' + str(disp[i]))
			outputFile.write('\n')
			

# get SS code from generated DSSP file
# input: file path
def getSS(input):
	secstrucs = []

	with open(input) as DSSP:
		for i, line in enumerate(DSSP):
			if i > 27: # skip header
				secstrucs.append(line[16]) # 1-letter ss code

	return secstrucs

# get disparities from previously generated file
# disparities: file path to disparities.txt
def getDisp(ID, disparities):
	with open(disparities) as disp:
		for i, line in enumerate(disp):
			if line[0:7] == ID:
				return line.split(',')[1:-1]

for i, ID in enumerate(IDs):
	if (ID != '2qe7G01') and (ID != '2qe7G02'):
		print(str(i) + ' ' + ID)
		generateCSV(ID)

# compiles data into 4-col table: ID,residueNum,ss,disp
# why ID? maybe compare superfamilies later
def SSversusDisp():
	with open('SSversusDisp.txt', 'w') as outputFile:
		for file in os.listdir('csv'):
			for i, line in enumerate(file):
				if i > 0: # skip header, which is ID
					outputFile.append(file.readlines()[0] + ',' + line[i])


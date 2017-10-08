import os
import subprocess

# clean sequence dataset -> ID only
# input file has ID and seq
def clean(i, o):
	clean = []

	with open(i) as f:
		for line in f.readlines():
			clean.append(line[0:7])

	with open(o, 'a') as f:
		for line in clean:
			f.write(line + '\n')

# use mkdssp executable
# NOTE: need dompbb folder, which is very large (download from CATH)
def generateDSSP(IDs):
	with open(IDs) as list:
		for ID in list.readlines(): # last char is \n
			subprocess.call('mkdssp -i dompdb/' + ID[:-1] + ' -o dssp/' + ID[:-1] + '.dssp', shell=True)

# put data into CSV to analyze
def generateCSV(dssp, csv):
	# resn,ss,disparity
	# sliding window is residue at resn + next 14

	ss = getSS(dssp)

	ID = dssp[5:-5]
	disp = getDisp(ID, 'disparities-n68.txt')

	with open(csv, 'a') as output:
		for i, item in enumerate(ss):
			if i in range(0, len(ss) - 17):
			# why -16? TODO: look into why # of disparity values is less than expected
				output.write(str(i) + ',' + item + ',' + disp[i] + '\n')
			else:
				output.write(str(i) + ',' + item + '\n')

# get SS code from generated DSSP file (file name is input)
def getSS(i):
	secstrucs = []

	with open(i) as f:
		for i, line in enumerate(f):
			if i > 27: # skip header
				secstrucs.append(line[16]) # 1-letter ss code

	return secstrucs

# get disparities from previously generated file
def getDisp(ID, disparities):
	with open(disparities) as f:
		for i, line in enumerate(f):
			if line[0:7] == ID:
				return str(next(enumerate(f))[1]).split()

# clean('final-dataset-n68.txt', 'IDs-n68.txt')
# generateDSSP('IDs-n68.txt')

# generate CSV for all 68
# NOTE: as of 10/7, it generates 65 and then list index out of range line 38
# current goal is to analyze preliminary data, so fix this later
with open('IDs-n68.txt') as IDs:
	for line in IDs:
		input = 'dssp/' + line[:-1] + '.dssp'
		output = 'csv/' + line[:-1] + '.csv'
		generateCSV(input, output)
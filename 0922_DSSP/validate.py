import os
import subprocess

# clean sequence dataset -> ID only
def clean(i, o):
	clean = []

	with open(i) as f:
		for line in f.readlines():
			clean.append(line[0:7])

	with open(o, 'a') as f:
		for line in clean:
			f.write(line + '\n')

# clean('final-dataset-n68.txt', 'IDs-n68.txt')

# use mkdssp executable
def generateDSSP(IDs):
	with open(IDs) as list:
		for ID in list.readlines(): # last char is \n
			subprocess.call('mkdssp -i dompdb/' + ID[:-1] + ' -o dssp/' + ID[:-1] + '.dssp', shell=True)

# generateDSSP('IDs-n68.txt')



# put data into CSV to analyze
def generateCSV(dssp, csv):
	# resn,ss,disparity
	# sliding window is residue at resn + next 14

	ss = getSS(dssp)

	ID = dssp[5:-5]
	disp = getDisp(ID, 'disparities-n68.txt')

	with open(csv, 'a') as output:
		for i, item in enumerate(ss):
			if disp[i]:
				output.write(str(i) + ',' + item + ',' + disp[i] + '\n')
			else:
				output.write(str(i) + ',' + item + '\n')

def getSS(i):
	secstrucs = []

	with open(i) as f:
		for i, line in enumerate(f):
			if i > 27: # skip header
				secstrucs.append(line[16]) # 1-letter ss code

	return secstrucs


def getDisp(ID, disparities):
	with open(disparities) as f:
		for i, line in enumerate(f):
			if line[0:7] == ID:
				return str(next(enumerate(f))[1]).split()

#print( getDisp('4fcyB02', 'disparities-n68.txt')[0] )

#getDisp('4fcyB02', 'disparities-n68.txt')

generateCSV('dssp/1qhdA02.dssp', 'csv/1qhdA02.csv')
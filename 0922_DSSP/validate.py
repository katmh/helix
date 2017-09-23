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

def generateDSSP(IDs):
	with open(IDs) as list:
		for ID in list.readlines(): # last char is \n
			subprocess.call('mkdssp -i dompdb/' + ID[:-1] + ' -o dssp/' + ID[:-1] + '.dssp', shell=True)

# generateDSSP('IDs-n68.txt')
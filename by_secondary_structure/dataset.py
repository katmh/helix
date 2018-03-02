from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
import math

def genbankID(i,o):
	with open(i) as f:
		with open(o, 'w') as out:
			for line in f:
				if line.startswith('GB:'):
					out.write(line)

#genbankID('ae_gb147_human_introns.flat', 'introns_gbIDs.txt')
#genbankID('ae_gb147_human_exons.flat', 'exons_gbIDs.txt')

def removeComp(i,o):
	with open(i) as f:
		with open(o, 'w') as out:
			for line in f:
				if line.split(' ')[-1].startswith('c') == False:
					out.write(line)

#removeComp('introns_gbIDs.txt', 'introns_gbIDs_noC.txt')
#removeComp('exons_gbIDs.txt', 'exons_gbIDs_noC.txt')

def thousand(inp,o,n):
	with open(inp) as f:
		with open(o, 'w') as out:
			for i, line in enumerate(f):
				if i % n == 0:
					out.write(line)

#thousand('introns_gbIDs_noC.txt', 'introns_gb_1000.txt', 51)
#thousand('exons_gbIDs_noC.txt', 'exons_gb_1000.txt', 34)

# then manually shave off the ends of each so there are 1000 entries

def download(i,o):
	Entrez.email = 'katherine_huang@student.uml.edu'

	with open(i) as f:
		with open(o, 'w') as out:
			for line in f:
				handle = Entrez.efetch(db='nucleotide', id=line[:-1], rettype='gb', retmode='text')
				record = SeqIO.read(handle, 'genbank')
				proteins = [f for f in record.features if f.type == 'Protein']

				for i, protein in enumerate(proteins):
					out.write(line[:-1] + ',') # write ID

					protein_left_limit = int(str(proteins[i].location).split(':')[0].lstrip('[').lstrip('<'))
					protein_right_limit = int(str(proteins[i].location).split(':')[1].rstrip(']').lstrip('>'))

					out.write(str(record.seq[protein_left_limit : protein_right_limit]) + '\n')

#download('exons_gb_1000.txt', 'exons_final.txt')
#download('introns_gb_1000.txt', 'introns_final.txt')

download('alpha_gene_IDs.list', 'alpha_proteins_sequences.txt')

# custom translate function to deal w/ non-multiples of 3
def trans(seq):
	if len(seq) % 3 == 0:
		return seq.translate()
	else:
		return seq[ : -(len(seq) % 3)].translate()

def disp(seq, size):
	table = {'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9, 'A': 1.8, 'G': -0.4, 'T': -0.7, 'W': -0.9, 'S': -0.8, 'Y': -1.3, 'P': -1.6, 'H': -3.2, 'E': -3.5, 'D': -3.5, 'N': -3.5, 'Q': -3.5, 'K': -3.9, 'R': -4.5, 'X': 0, '*': 0, 'J': 4.5, 'Z': -3.5}

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

# exons are denoted 1 (positive example), introns 0 (negative)

def makeDataset():
	with open('dataset_disp.csv', 'w') as out:
		out.write('seq,ie,gc,nstart,maxnstop,maxdisp\n')
		with open('exons_final.txt') as f:
			for line in f:
				# I/E label
				out.write(line[:-1] + ',1,')

				# GC
				my_seq = Seq(line[:-1], IUPAC.unambiguous_dna)
				out.write(str(GC(my_seq)) + ',')

				# N_ATG (start codon)
				out.write(str(my_seq.count_overlap('ATG')) + ',')

				# max(N_TAA, N_TAG, N_TGA) (stop codons)
				out.write(str( max( my_seq.count_overlap('TAA'), my_seq.count_overlap('TAG'), my_seq.count_overlap('TGA') ) ) + ',')

				# max disparity in translation of default frame
				#out.write( str(max(disp(trans(my_seq), 15))) + '\n')
		with open('introns_final.txt') as f:
			for line in f:
				# I/E label
				out.write(line[:-1] + ',0,')

				# GC
				my_seq = Seq(line, IUPAC.unambiguous_dna)
				out.write(str(GC(my_seq)) + ',')

				# N_ATG (start codon)
				out.write(str(my_seq.count_overlap('ATG')) + ',')

				# max(N_TAA, N_TAG, N_TGA) (stop codons)
				out.write(str( max( my_seq.count_overlap('TAA'), my_seq.count_overlap('TAG'), my_seq.count_overlap('TGA') ) ) + ',')

				# max disparity in translation of default frame
				#out.write( str(max(disp(trans(my_seq), 15))) + '\n')

#makeDataset()

def transFile():
	with open('translated.txt', 'w') as out:
		with open('dataset.csv') as f:
			for i, line in enumerate(f):
				if i != 0:
					my_seq = Seq(line.split(',')[0], IUPAC.ambiguous_dna)
					out.write(str(trans(my_seq)) + '\n')

#transFile()

def dispFile():
	with open('disparities.txt', 'w') as out:
		with open('translated.txt') as trans:
			for line in trans:
				for val in disp(line[:-1], 15):
					out.write(str(val) + ',')
				out.write('\n')

#dispFile()

def maxDispFile():
	with open('disparities.txt') as f:
		with open('maxdisparities.txt', 'w') as out:
			for i, line in enumerate(f):
				if len( line.split(',')[:-1] ) == 0: # if not blank
					out.write('NA\n')
					print('zero at ' + str(i))
				else:
					out.write( str(max(line.split(',')[:-1])) + '\n' )

#maxDispFile()

'''
19 NAs in disp:
zero at 23
zero at 75
zero at 78
zero at 133
zero at 175
zero at 304
zero at 351
zero at 493
zero at 593
zero at 778
zero at 811
zero at 879
zero at 975
zero at 1004
zero at 1007
zero at 1261
zero at 1381
zero at 1383
zero at 1436

TODO: remove sequences containing mostly Ns
'''
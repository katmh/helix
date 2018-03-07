from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
import math
import statistics
import itertools

def getPDBIDs(i, o):
	print('extracting IDs from PDB info file ' + i + ' into ' + o)
	with open(i) as f:
		with open(o, 'w') as out:
			for i, line in enumerate(f):
				if i != 0:
					out.write(line.split(',')[0].strip('"') + '\n')

#getPDBIDs('alpha_proteins_PDB_info.csv', 'alpha_proteins_IDs.list')

#getPDBIDs('beta_proteins_PDB_info.csv', 'beta_proteins_IDs.list')
#getPDBIDs('ab_proteins_PDB_info.csv', 'ab_proteins_IDs.list')
#getPDBIDs('fewSS_proteins_PDB_info.csv', 'fewSS_proteins_IDs.list')
#getPDBIDs('beta_test.csv', 'beta_test_protein_IDs.list') # test with 30% seq identityn

# UniProt IDs mapped to multiple RefSeq IDs; I only want 1 RefSeq ID per UniProt ID
def removeDuplicateIDs(i, o):
	with open(i) as f:
		with open(o, 'w') as out:
			uniprotIDs = list()
			for line in f:
				if line.split('	')[0] not in uniprotIDs:
					uniprotIDs.append(line.split('	')[0])
					out.write(line.split('	')[1])

#removeDuplicateIDs('alpha_uniprot_and_refseq_IDs.list', 'alpha_refseq_IDs.list')

def download(i,o):
	print('downloading ' + i + ' into ' + o)
	Entrez.email = 'katherine_huang@student.uml.edu'

	with open(i) as f:
		with open(o, 'w') as out:
			for n, line in enumerate(f):
				print(str(i) + ': line ' + str(n))
				handle = Entrez.efetch(db='nucleotide', id=line[:-1], rettype='gb', retmode='text')
				record = SeqIO.read(handle, 'genbank')

				CDSs = [f for f in record.features if f.type == 'CDS']

				for j, CDS in enumerate(CDSs):
					# out.write(line[:-1] + ',') # write ID
					if j == 0:
						CDS_left_limit = str(CDS.location).split(':')[0].lstrip('[')
						if len(CDS_left_limit.split('[')) > 1:
							CDS_left_limit = int(CDS_left_limit.split('[')[1])
						if len(str(CDS_left_limit).split('<')) > 1:
							CDS_left_limit = int(CDS_left_limit.split('<')[1])
						else:
							CDS_left_limit = int(CDS_left_limit)
						CDS_right_limit = int(str(CDS.location).split(':')[1].split(']')[0])
						
						out.write(str(record.seq[CDS_left_limit : CDS_right_limit]) + '\n')

#download('alpha_refseq_IDs.list', 'alpha_coding.csv')
#download('beta_refseq_IDs.list', 'beta_coding.csv')
#download('ab_refseq_IDs.list', 'ab_coding.csv')
#download('fewSS_refseq_IDs.list', 'fewSS_coding.csv')

''' FILTERING '''
# filter out sequences that are >50000bp or <45 bp
# filter out sequences with 'N' bases
def filter(i, o):
	with open(i) as f:
		with open(o, 'w') as out:
			print('filtering ' + i + ' for length (must be <50000bp oand >45bp) and repeats (remove sequences with "N"')
			for line in f:
				# check length, skip newline char
				if len(line[:-1]) < 50000 and len(line[:-1]) > 45 and 'N' not in line[:-1]:
					out.write(line)

#filter('alpha_coding.csv', 'alpha_coding_filtered.csv')
#filter('beta_coding.csv', 'beta_coding_filtered.csv')
#filter('ab_coding.csv', 'ab_coding_filtered.csv')
#filter('fewSS_coding.csv', 'fewSS_coding_filtered.csv')
#filter('introns_final_8179.txt', 'introns_final_filtered.csv')

''' FEATURES '''

# custom translate function; need translations for disparity calculation
# if not multiple of 3, remove the remainder
def trans(seq):
	if len(seq) % 3 == 0:
		return seq.translate()
	else:
		return seq[ : -(len(seq) % 3)].translate()

def disp(seq):
	seq = trans(seq) # must be Seq object

	table = {'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9, 'A': 1.8, 'G': -0.4, 'T': -0.7, 'W': -0.9, 'S': -0.8, 'Y': -1.3, 'P': -1.6, 'H': -3.2, 'E': -3.5, 'D': -3.5, 'N': -3.5, 'Q': -3.5, 'K': -3.9, 'R': -4.5, 'X': 0, '*': 0, 'J': 4.5, 'Z': -3.5}

	disparities = []

	# iterate thru seq and get starting residues of each window
	for i, start in enumerate( seq[0 : len(seq) - (15 - 1)] ):

		# variables to be populated
		hydro = []
		sumX = 0
		sumY = 0

		# from start residue, iterate thru next 15-1 residues
		for n, residue in enumerate( seq[i : i + 15] ):

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
		disparities.append(disp / 15)

	return statistics.median(disparities)

def ORF(seq): # todo: also add ORF coverage
	ORF_lengths = []
	for i, base in enumerate(str(seq)): # input was Seq object
		if i < len(seq) - 2:
			if seq[i] + seq[i+1] + seq[i+2] == 'ATG':
				stop_indices = []
				stop_codons = ['TAA', 'TAG', 'TGA']
				for codon in stop_codons:
					stop_indices.append(seq.find(codon, i+3))
				for index in stop_indices:
					if (i % 3 == index % 3) and (index != -1):
						ORF_lengths.append((index+3)-i)
	if len(ORF_lengths) > 0:
		return str(max(ORF_lengths)) + ',' + str(max(ORF_lengths)/len(seq))
	else:
		return '0,0'
'''
def hex(seq):
	bases = ['A', 'T', 'C', 'G']
	hexamers_list = list(itertools.product(bases, repeat=6))

'''
def extractFeature(i, o, feature):
	print('extracting ' + str(feature) + ' for ' + i)
	with open(i) as f:
		with open(o, 'w') as out:
			for i, line in enumerate(f):
				print('line ' + str(i))
				if i != 0:
					out.write(str(feature(Seq(line.split(',')[1].strip('"')))) + '\n')

#extractFeature('data_full.csv', 'data_full_disp.csv', disp)
#extractFeature('data_full.csv', 'data_full_GC.csv', GC)
#extractFeature('data_full.csv', 'data_full_ORF.csv', ORF)
#extractFeature('data.csv', 'data_disp.csv', disp)
#extractFeature('beta_coding.csv', 'beta_coding_disp.csv', disp)
#extractFeature('ab_coding.csv', 'ab_coding_disp.csv', disp)
#extractFeature('fewSS_coding.csv', 'fewSS_coding_disp.csv', disp)
#extractFeature('data.csv', 'data_GC.csv', GC)
#extractFeature('beta_coding.csv', 'beta_coding_GC.csv', GC)
#extractFeature('ab_coding.csv', 'ab_coding_GC.csv', GC)
#extractFeature('fewSS_coding.csv', 'fewSS_coding_GC.csv', GC)

#extractFeature('data.csv', 'data_ORF.csv', maxORF)

with open('data_full_ORF.csv') as f:
	with open('data_full_ORFcover.csv', 'w') as out:
		for line in f:
			out.write(line.split(',')[1])
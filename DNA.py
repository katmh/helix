# input should be template 3' -> 5', output mRNA will be 5' -> 3', which is the direction of translation
# parameters: input file of DNA w/ header
def transcribe(i):
	DNAseq = ''

	with open(i) as DNA:
		for i, line in enumerate(DNA):
			# skip header
			if i != 0:
				DNAseq += line[:-1] # don't include newline char

	dnaRna = {'C': 'G', 'G': 'C', 'A': 'U', 'T': 'A'}

	RNA = ''

	for char in DNAseq:
		RNA += dnaRna[char]

	return RNA

#print(translate('TCGTA'))

# parameter: input transcript (nucleotide FASTA w/ header)
def translate(transcript):
	# one-letter amino acids
	codons = {
		# U
		'UUU': 'F',
		'UUC': 'F',
		'UUA': 'L',
		'UUG': 'L',

		'UCU': 'S',
		'UCC': 'S',
		'UCA': 'S',
		'UCG': 'S',

		'UAU': 'Y',
		'UAC': 'Y',
		'UAA': '-', # stop
		'UAG': '-', # stop

		'UGU': 'C',
		'UGC': 'C',
		'UGA': '-', # stop
		'UGG': 'W',

		# C
		'CUU': 'L',
		'CUC': 'L',
		'CUA': 'L',
		'CUG': 'L',

		'CCU': 'P',
		'CCC': 'P',
		'CCA': 'P',
		'CCG': 'P',

		'CAU': 'H',
		'CAC': 'H',
		'CAA': 'Q',
		'CAG': 'Q',

		'CGU': 'R',
		'CGC': 'R',
		'CGA': 'R',
		'CGG': 'R',

		# A
		'AUU': 'I',
		'AUC': 'I',
		'AUA': 'I',
		'AUG': 'M',

		'ACU': 'T',
		'ACC': 'T',
		'ACA': 'T',
		'ACG': 'T',

		'AAU': 'N',
		'AAC': 'N',
		'AAA': 'K',
		'AAG': 'K',

		'AGU': 'S',
		'AGC': 'S',
		'AGA': 'R',
		'AGG': 'R',

		# G
		'GUU': 'V',
		'GUC': 'V',
		'GUA': 'V',
		'GUG': 'V',

		'GCU': 'A',
		'GCC': 'A',
		'GCA': 'A',
		'GCG': 'A',

		'GAU': 'D',
		'GAC': 'D',
		'GAA': 'E',
		'GAG': 'E',

		'GGU': 'G',
		'GGC': 'G',
		'GGA': 'G',
		'GGG': 'G'
	}

	# amino acid sequence to be returned
	aaSeq = []

	# split RNA into codons and find corresponding amino acid
	for i, base in enumerate(transcript):
		if i % 3 == 0 and i < len(transcript) - 2: # let base at i be the start of codon; cannot be last 2 bases
			print(i)
			aaSeq.append( codons[base + transcript[i+1] + transcript[i+2]] )

	return aaSeq

print(translate(transcribe('DPYD-transcript-variant-1.txt')))
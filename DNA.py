# input should be template 3' -> 5', output mRNA will be 5' -> 3', which is the direction of translation
def transcribe(DNA):
	dnaRna = {'C': 'G', 'G': 'C', 'A': 'U', 'T': 'A'}
	RNA = ''
	for char in DNA:
		RNA += dnaRna[char]
	return RNA

# input is mRNA 5' -> 3'
def translate(transcript):
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
	aaSeq = ''
	for i, base in enumerate(transcript):
		if i % 3 == 0 and i < len(transcript) - 2: # let base at i be the start of codon; cannot be last 2 bases
			aaSeq += codons[base + transcript[i+1] + transcript[i+2]] # match codon
	return aaSeq

#print(translate(transcribe('DPYD-transcript-variant-1.txt')))

with open('NONCODE2016_human_clean_every20.txt') as f:
	with open('seq_NONCODE.csv', 'w') as out:
		lines = f.readlines()
		for i, line in enumerate(lines):
			if line.startswith('>'):
				out.write(line[1:-1] + ',')
				print('line ' + str(i) + ' written to output')

				out.write( translate( transcribe( lines[i+1].rstrip('\n').upper() ) ) + '\n' )
				print('line ' + str(i+1) + ' transcribed, translated, and written')
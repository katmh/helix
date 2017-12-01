# input should be template 3' -> 5', output mRNA will be 5' -> 3', which is the direction of translation
def transcribe(DNA):
	dnaRna = {'C': 'G', 'G': 'C', 'A': 'U', 'T': 'A'}

	RNA = ''

	for char in DNA:
		RNA += dnaRna[char]

	return RNA

#print(translate('TCGTA'))

def translate(RNA):
	# one-letter amino acids
	codons = {}

	# split RNA into 3s

	return seq
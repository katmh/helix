# parameter: i is FASTA file
def cleanFile(i, o):
	with open(i) as f:
		with open(o, 'w') as out:
			for i, line in enumerate(f):
				if line[0] == '>':
					out.write('\n' + line) # including \n
					print('wrote line ' + str(i) + ' with newline in front')
				elif line[0] != '>':
					out.write(line.rstrip('\n'))
					print('wrote line ' + str(i) + ' without newline')

#cleanFile('NONCODE2016_human.fa', 'NONCODE2016_human_clean.txt')

def every20(i, o):
	with open(i) as f:
		with open(o, 'w') as out:
			lines = f.readlines()
			for i in range(0, len(lines)):
				if i % 40 == 0: # mod 40 bc each entity takes 2 lines
					out.write(lines[i] + lines[i+1])

#every20('NONCODE2016_human_clean.txt', 'NONCODE2016_human_clean_every20.txt')


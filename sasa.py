# verify that all PDB IDs in CSV exist in the SASA masterlist
def verifyExist(csv, sasa):
	with open(csv) as CSV:
		with open(sasa) as SASA:
			sasaLines = SASA.readlines()
			sasaIDs = []
			for i, entry in enumerate(sasaLines):
				if i != 0:
					sasaIDs.append(entry[2:6].upper())

			numFound = 0

			for i, line in enumerate(CSV):
				if i != 0: # skip header
					if line.split(',')[0].strip('"') not in sasaIDs:
						print(line.split(',')[0].strip('"') + ' not found')
					if line.split(',')[0].strip('"') in sasaIDs:
						numFound += 1

			print(str(numFound) + ' chains found')

#verifyExist('seq_mainly_beta.csv', 'pdb_sasa.txt')

def phobicPhilic(i, o):
	# pos is hydrophobic, neg is hydrophilic
	table = {'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9, 'A': 1.8, 'G': -0.4, 'T': -0.7, 'W': -0.9, 'S': -0.8, 'Y': -1.3, 'P': -1.6, 'H': -3.2, 'E': -3.5, 'D': -3.5, 'N': -3.5, 'Q': -3.5, 'K': -3.9, 'R': -4.5, 'X': 0, 'U': 2.5, 'O': -3.9, 'B': -3.5, '-': 0}

	with open(i) as CSV:
		with open(o, 'w') as out:
			out.write('ID,percentPhobic,percentPhilic\n')

			for i, line in enumerate(CSV):
				if i != 0:
					# write ID
					out.write(line.split(',')[0].strip('"') + ',')

					# write percent hydrophobic
					numPhobic = 0
					seq = line.split(',')[3].split('#')[0].lstrip('"')

					for residue in seq:
						if table[residue] > 0:
							numPhobic += 1

					out.write( str( float( (numPhobic)/(len(seq)) ) ) + ',' )

					# write percent hydrohpilic
					out.write( str( float( (len(seq)-numPhobic)/(len(seq)) ) ) + '\n' )

#phobicPhilic('seq_mainly_beta.csv', 'pp_mainly_beta.csv')
#phobicPhilic('seq_mainly_alpha.csv', 'pp_mainly_alpha.csv')
phobicPhilic('seq_alpha_beta.csv', 'pp_ab.csv')

def SASA(i, sasa, o):
	with open(i) as CSV:
		with open(sasa) as SASA:
			sasaLines = SASA.readlines()[2:]

			with open(o, 'w') as out:
				# write CSV header
				out.write('ID,normSASA\n')

				for i, line in enumerate(CSV):
					if i != 0:
						# write ID and comma
						out.write(line.split(',')[0].strip('"') + ',')

						# normalize SASA: area / # atoms (bc obviously a larger protein might have more SASA)
						normSASA = 0.0
						for entry in sasaLines:
							if line.split(',')[0].strip('"') == entry[2:6].upper():
								normSASA += float(entry[13:23].strip(' ')) / float(entry[6:13].strip(' '))

								# write normalized SASA
								out.write(str(normSASA) + '\n')
								print(str(i) + ': ' + str(normSASA) + ' written to file')

#SASA('seq_mainly_beta.csv', 'pdb_sasa.txt', 'sasa_mainly_beta.csv')
#SASA('seq_mainly_alpha.csv', 'pdb_sasa.txt', 'sasa_mainly_alpha.csv')
SASA('seq_alpha_beta.csv', 'pdb_sasa.txt', 'sasa_ab.csv')
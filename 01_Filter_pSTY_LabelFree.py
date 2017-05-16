import sys; import urllib2; import numpy; import os
def readTab(infile): # read in txt file
	with open(infile, 'r') as input_file: 
	# read in tab-delim text
		output = []
		for input_line in input_file:
			input_line = input_line.strip()
			temp = input_line.split('\t')
			output.append(temp)
	return output
def get_info(uniprot_accession_in): #get aa lengths and gene name
	error = open('error proteins.txt', 'a+')
	i=0
	while i==0:
		try:
			data = urllib2.urlopen("http://www.uniprot.org/uniprot/" + uniprot_accession_in + ".fasta")
			break
		except urllib2.HTTPError, err:
			i = i + 1
			if i == 50:
				sys.exit("More than 50 errors. Check your file or try again later.")
			if err.code == 404:
				error.write(uniprot_accession_in + '\t' + "Invalid URL. Check protein" + '\n')
				seqlength = 'NA'
				genename = 'NA'
				return genename
			elif err.code == 302:
				sys.exit("Request timed out. Check connection and try again.")
			else:
				sys.exit("Uniprot had some other error")
	lines = data.readlines()
	header = lines[0]
	lst = header.split('|')
	lst2 = lst[2].split(' ')
	swissprot = lst2[0]
	uniprot_acc = lst[1]
	if lines == []:
		error.write(uniprot_accession_in + '\t' + "Blank Fasta" + '\n')
		error.close
		uniprot_acc = 'NA'
		genename = 'NA'
		return genename
	if lines != []:
		seqlength = 0
		header = lines[0]
		if 'GN=' in header:
			lst = header.split('GN=')
			lst2 = lst[1].split(' ')
			genename = lst2[0]
			error.close
			return genename
		if 'GN=' not in header:
			genename = 'NA'
			error.close
			return genename
def FilterMaxquant(pSTY_file, PEP_cutoff, MZ_error_cutoff):
	data = readTab(pSTY_file)
	PEP_index = data[0].index("PEP")
	Int_index = data[0].index("Intensity")
	MZerror_index = data[0].index("Mass Error [ppm]")
	samples = []
	for i in data[0]:
		if i.startswith("Intensity "):
			if i != "Intensity":
				samples.append(i)
		aa_index = data[0].index("Amino Acid")
		Position_index = data[0].index("Position")
		Protein_index = data[0].index("Protein")
		modSeq_index = data[0].index("Modified Sequence")
		filtered = [["Source","ID","Accession","Gene","PEP","Mass Error [ppm]"]]
	for i in samples:
		filtered[0].append(i)
	for i in data[1:]:
		temp = ["Label Free"]
		#print i[PEP_index]
		if float(i[PEP_index]) < float(PEP_cutoff):
			if abs(float(i[MZerror_index])) < float(MZ_error_cutoff):
				if float(i[Int_index]) != 0:
					ID = i[Protein_index] + i[modSeq_index] + i[aa_index] + i[Position_index]
					temp.append(ID);temp.append(i[Protein_index]); temp.append(get_info(i[Protein_index]))
					temp.append(i[PEP_index]);temp.append(i[MZerror_index])
					for j in samples:
			   			temp.append(i[data[0].index(j)])
					filtered.append(temp)
	return filtered

def writeTable(DF):
	with open("pY_output_automated.txt",'w') as x:
		for i in DF:
			x.write('\t'.join(i));x.write("\n")

writeTable(FilterMaxquant(sys.argv[1], sys.argv[2], sys.argv[3]))
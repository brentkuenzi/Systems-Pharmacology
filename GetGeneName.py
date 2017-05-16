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
def writeTable(DF):
	with open("tmpDF2.txt",'w') as x:
		for i in DF:
			x.write('\t'.join(i));x.write("\n")

data = readTab(sys.argv[1])
ID = data[0].index(sys.argv[2])
data[0].append(sys.argv[3])
for i in data[1:]:
	i.append(get_info(i[ID]))
writeTable(data)
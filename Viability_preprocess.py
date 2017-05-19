import sys
import numpy
def readTab(infile): # read in txt file
    with open(infile, 'r') as input_file:
    # read in tab-delim text
        output = []
        for input_line in input_file:
            input_line = input_line.strip()
            temp = input_line.split('\t')
            output.append(temp)
    return output
def backgroundSubtract(table,row): # calculate and remove background
	sum = 0
	output=[]
	for i in table[row-1]:
		sum += float(i)
	background = sum/len(table[row-1])
	for i in table:
		temp = []
		for j in i:
			temp.append(str(float(j)-background))
		output.append(temp)
	return output
def normalizeViability(table, row=False,col=False): # normalize counts to vehicle control row or column
	if row != False:
		sum = 0
		output=[]
		for i in table[row-1]:
			sum += float(i)
		vehicle = sum/len(table[row-1])
		for i in table:
			temp = []
			for j in i:
				temp.append(str(100*(float(j)/vehicle)))
			output.append(temp)
		return output
	if col != False:
		sum = 0
		cnt = 0
		output=[]
		for i in table:
			sum += float(i[col-1])
			cnt+=1
		vehicle = sum/cnt
		for i in table:
			temp = []
			for j in i:
				temp.append(str(100*(float(j)/vehicle)))
			output.append(temp)
		return output
def normalizeViabilityLowest(table,dose,by_drug = "None",desc=False):
	normalizer = 1
	if by_drug != "None":
		for i in dose[1:]:
			if by_drug == i[4]:
				temp = []
				temp = i[0] + "-" + i[1]
				positions = temp.split("-")
		temp = []
		for i in range(0,len(positions)):
			positions[i] = int(positions[i])-1
		for j in range(positions[0],positions[1]+1):
			temp.append(float(table[j][positions[2]]))
		normalizer = numpy.mean(temp,axis=0)
	for i in dose[1:]: # row start, row end, col start, col end
		temp = []
		temp = i[0] + "-" + i[1]
		positions = temp.split("-")
		for i in range(0,len(positions)):
			positions[i] = int(positions[i])-1
		if desc == False:
			temp = []
			for j in range(positions[0],positions[1]+1):
				temp.append(float(table[j][positions[2]]))
			if by_drug == "None":
				normalizer = numpy.mean(temp,axis=0)
			for j in range(positions[0],positions[1]+1):
				for k in range(positions[2],positions[3]+1):
					table[j][k] = str(100*(float(table[j][k])/normalizer))
		else:
			temp = []
			for j in range(positions[0],positions[1]+1):
				temp.append(float(table[j][positions[3]]))
			normalizer = numpy.mean(temp)
			for j in range(positions[0],positions[1]+1):
				for k in range(positions[2],positions[3]+1):
					table[j][k] = str(100*(float(table[j][k])/normalizer))

	return table
def doseResponse(data_file,dose_file,background_row,by_drug = "None", desc=False,lowest0 = True,normalize="lowest"):
	# data_file =  file you want to upload
	# dose_file = user generated info file about plating
	# background_row = row to use for background correction
	# format = M5 spectramax output or simple matrix format
	# desc = is the data set high to low conc?
	# lowest0 = is the lowest dose vehicle?
	# normalize options = "lowest",["row",#],["col",#]
		# normalize by lowest concentration for each curve or normalize
		# everything by a certain column?

	data = readTab(data_file)
	dose = readTab(dose_file)
	doses = {}
	location = {}
	drugs = []
	for i in dose[1:]: # get doses for each drug and store as dictionary
		temp = []
		temp = i[0] + "-" + i[1]
		positions = temp.split("-")
		for j in range(0,len(positions)):
			positions[j] = int(positions[j])-1
		location[i[4]] = positions
		drugs.append(i[4])
		cnt=0
		temp = []
		if "HighConc" in dose[0]:
			for j in range(positions[2],positions[3]):
				if cnt == 0:
					temp.append(float(i[2]))
				else:
					temp.append(temp[cnt-1]/(int(i[3])))
				cnt+=1
			doses[i[4]] = temp
		else:
			doses[i[4]] = i[2].split(",")
	data_bg = backgroundSubtract(data, background_row) # subtract background
	fullData = {}
	averages = {}
	sd = {}
	if normalize != "lowest": # normalization method
		if normalize[0] == "col":
			data_nm = normalizeViability(data_bg, col=normalize[1])
		if normalize[0] == "row":
			data_nm = normalizeViability(data_bg, row=normalize[1])
	else:
		data_nm = normalizeViabilityLowest(data_bg, dose,by_drug = by_drug, desc=desc)


	with open("output.txt","w") as x: # write file to import into R
		x.write("\t".join(["Dose","Response","Drug"]))
		x.write("\n")
		for i in drugs:
			for j in range(location[i][0],location[i][1]+1):
					cnt = 0
					for k in range(location[i][2],location[i][3]+1):
						if desc == False:
							if lowest0 == True:
								if cnt == 0:
									x.write("0")
								else:
									x.write(str(doses[i][-1*(cnt)]))
						else:
							if lowest0 == True:
								if cnt == len(doses[i]):
									x.write("0")
							else:
								x.write(str(doses[i][cnt]))
						x.write("\t"); x.write(str(data_nm[j][k])); x.write("\t")
						x.write(i); x.write("\n")
						cnt+=1
	# organize into 3 column table for each drug (dose, response, drug)
	# play with drc to get multiple lines on same plot
	# will then add points w/ error bars for each drug
doseResponse(sys.argv[1], sys.argv[2], 2, by_drug = sys.argv[3], desc = False, lowest0 = True, normalize="lowest")

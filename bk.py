import sys; import os
def convert_units(value,input_format="M",output_format="uM"):
	levels = {"M":float(1),"mM":float(1000), "uM":float(1000000), "nM":float(1000000000),
	"L":float(1),"mL":float(1000),"uL":float(1000000),"nL":float(1000000000),
	"liter":float(1),"milliliter":float(1000),"microliter":float(1000000),"nanoliter":float(1000000000),
	"molar":float(1),"millimolar":float(1000),"micromolar":float(1000000),"nanomolar":float(1000000000),
	"g":float(1),"mg":float(1000),"ug":float(1000000),"ng":float(1000000000),
	"grams":float(1),"milligrams":float(1000),"micrograms":float(1000000),"nanograms":float(1000000000)}
	return float(value / levels[input_format]) * levels[output_format]
def calc_molarity(weight,weight_units,mw,volume,volume_units):
	return float((convert_units(weight,weight_units,"g")/mw)/convert_units(volume,volume_units,"L"))
def calc_dilution(C1,C1_units,V1,V1_units,C2,C2_units,V2,V2_units,output_format = "value"):
	flag = 0
	for i in [C1,C2,V1,V2]:
		if flag > 1:
			sys.exit("ERROR: More than 1 argument is missing!")
		if i == None:
			flag+=1
	if output_format not in ["value","string"]:
		sys.exit("ERROR: Input valid output_format")
	if C1 == None:
		value = (float(convert_units(C2,C2_units,"M")) * convert_units(V2,V2_units,"L")) / convert_units(V1,V1_units,"L")
		if output_format == "value":
			return convert_units(value,"M",C1_units)
		if output_format == "string":
			return str(convert_units(value,"M",C1_units)) + C1_units
	if C2 == None:
		value = (float(convert_units(C1,C1_units,"M")) * convert_units(V1,V1_units,"L")) / convert_units(V2,V2_units,"L")
		if output_format == "value":
			return convert_units(value,"M",C2_units)
		if output_format == "string":
			return str(convert_units(value,"M",C2_units)) + C2_units
	if V1 == None:
		value = (float(convert_units(C2,C2_units,"M")) * convert_units(V2,V2_units,"L")) / convert_units(C1,C1_units,"M")
		if output_format == "value":
			return convert_units(value,"L",V1_units)
		if output_format == "string":
			return str(convert_units(value,"L",V1_units)) + V1_units
	if V2 == None:
		value = (float(convert_units(C1,C1_units,"M")) * convert_units(V1,V1_units,"L")) / convert_units(C2,C2_units,"M")
		if output_format == "value":
			return convert_units(value,"L",V2_units)
		if output_format == "string":
			return str(convert_units(value,"L",V2_units)) + V2_units

#print "Conversion test..." + "Passed!" if convert_units(1,"M","mM") == 1000 else "Failed!"
#print "Dilution test..." + "Passed!" if calc_dilution(C1 = 10,C1_units = "mM",V1 = None,V1_units="uL",C2 = 1.5,C2_units="uM",V2=3,V2_units="mL",output_format="string") == "0.45uL" else "Failed!"
#print "Molarity test..." + "Passed!" if calc_molarity(531.45,"g",531.45,1,"L")==1 else "Failed"
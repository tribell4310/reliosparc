"""

Tristan Bell
Chao Lab
Massachusetts General Hospital

This script takes a cryosparc .cs file for particles or particle passthrough and exports the particle coordinates to a
group of star files that can mimic a relion autopicking job.

Script assumes that your motion-corrected micrographs for extraction are in mrc format.


"""

import sys
import numpy as np
import os
from os import listdir
from os.path import isfile, join


def main(inCs):
	# Check for the Raw_data subdirectory and delete any contents
	if os.path.isdir("./Raw_data") == False:
		os.mkdir("./Raw_data")
	else:
		print("\nDetected a Raw_data folder in this directory.")
	onlyfiles = [f for f in listdir("./Raw_data/") if isfile(join("./Raw_data", f))]
	if len(onlyfiles) > 0:
		response = input("Files detected in the Raw_data folder. Okay to overwrite them? [y/n] ")
		if response in ["Y", "y"]:
			print("Okay.  Deleting files in Raw_Data...")
			for item in onlyfiles:
				os.remove(join("./Raw_data/", item))
		else:
			print("Exiting.")
			exit()

	# Read in the star file as a np array
	print("Reading the cryosparc file...")
	f = np.load(inCs)

	# Get the starting index for the particle location info
	print("Procesing micrograph names...")
	start_index = infer_index(f)

	# Sample the micrograph names to find the invariant portions
	if len(f) < 100000:
		end_ind = len(f)
	else:
		end_ind = 100000

	mcg_names = []
	for i in range(0, end_ind-1):
		mcg_names.append(f[i][start_index].decode("utf-8"))
	
	# Work from the first and last non-constant indeces to identify the micrograph names
	constant_container = get_constant_matrix(mcg_names)
	first_var = 0
	for i in range(0, len(constant_container)):
		if constant_container[i] == False:
			if first_var == 0:
				first_var = i
			last_var = i
	mcg_start = mcg_names[0].find("_", first_var) + 1
	mcg_suffix = mcg_find_suffix(mcg_names, mcg_start)
	
	# Load the particles into a dictionary
	print("Transforming particle coordinates...")
	coord_dict = {}
	counter = 1
	for item in f:
		# Get micrograph name and confirm in dictionary
		filename = item[start_index].decode("utf-8")
		end_index = filename.find(mcg_suffix)
		mcg_name = filename[mcg_start:end_index]+".mrc"
		if counter % 100000 == 0:
			print(str(counter)+" / "+str(len(f)))
		counter += 1
		if mcg_name not in coord_dict:
			coord_dict[mcg_name] = {"x":[], "y":[]}

		# Calculate transformed x and y coords
		h = item[start_index+1][0]
		l = item[start_index+1][1]
		x_frac = item[start_index+2]
		y_frac = item[start_index+3]
		x_coord = round((l*x_frac), 0)
		y_coord = round(h-(h*y_frac), 0)

		# Add to dictionary
		coord_dict[mcg_name]["x"].append(str(round(x_coord, 6)))
		coord_dict[mcg_name]["y"].append(str(round(y_coord, 6)))

	# Write out the coords as "autopicking" star files
	print("Writing star files...")
	counter = 1
	for item in coord_dict.keys():
		if counter % 1000 == 0:
			print(str(counter)+" / "+str(len(coord_dict.keys())))
		counter += 1
		g = open("Raw_data/"+no_ext(item)+"_autopick.star", "w", newline="")
		#g.write("\n# version 30001\n\ndata_particle\n\nloop_ \n_rlnCoordinateX #1 \n_rlnCoordinateY #2 \n_rlnAutopickFigureOfMerit #3 \n_rlnClassNumber #4 \n_rlnAnglePsi #5 \n")
		g.write("\n# version 30001\n\ndata_\n\nloop_ \n_rlnCoordinateX #1 \n_rlnCoordinateY #2 \n_rlnAutopickFigureOfMerit #3 \n_rlnClassNumber #4 \n_rlnAnglePsi #5 \n")
		for i in range(0, len(coord_dict[item]["x"])):
			g.write(line_writer(coord_dict[item]["x"][i], coord_dict[item]["y"][i]))
		g.write("\n")
		g.close()
	print("Done.")


def mcg_find_suffix(full_list, start_ind):
	# Take the substring from the start of mcg name to end of the full cryosparc name
	names_list = []
	for i in range(0, len(full_list)):
		names_list.append(full_list[i][start_ind:])
	
	# Leftpad all names to the same length
	max_len = 0
	for i in range(0, len(names_list)):
		if len(names_list[i]) > max_len:
			max_len = len(names_list[i])
	same_len_names = []
	for i in range(0, len(names_list)):
		temp_str = names_list[i]
		if len(names_list[i]) < max_len:
			while len(temp_str) < max_len:
				temp_str = " "+temp_str
		same_len_names.append(temp_str)

	# Make an invariance matrix
	const_matrix = get_constant_matrix(same_len_names)

	# Working backwards, find the first True-False transition
	for i in range(1, len(const_matrix)):
		if const_matrix[-i] == True:
			if const_matrix[-(i+1)] == False:
				end_index = -i
				break

	return same_len_names[0][end_index:]


def get_constant_matrix(mcg_names):
	constant_container = []
	for i in range(0, len(mcg_names[0])):
		is_constant = True
		for j in range(1, len(mcg_names)):
			try:
				if len(mcg_names[j][i]) == len(mcg_names[0][i]):
					if len(mcg_names[j-1][i]) == len(mcg_names[0][i]):
						if mcg_names[j][i] != mcg_names[j-1][i]:
							is_constant = False
							break
			except:
				pass
		constant_container.append(is_constant)
	
	return constant_container

def infer_index(np_array):
	# Defined pattern is binary string, list of two ints >1000, float <= 1, float <=1
	for i in range(0, len(np_array[1])):
		try:
			np_array[1][i].decode("utf-8")
			try:
				a = len(np_array[1][i+1])
				if np_array[1][i+2] < 1:
					if np_array[1][i+3] < 1:
						startInd = i
						break
			except:
				pass
		except:
			pass

	return startInd


def line_writer(x, y):
	# Process x and y
	padded_x = leftpad(x, 12)
	padded_y = leftpad(y, 12)
	remainder = "     0.080000            0     0.000000 \n"
	return (padded_x + " " + padded_y + remainder)


def leftpad(inStr, final_len):
	while len(inStr) < final_len:
		inStr = " "+inStr
	return inStr


def no_ext(inStr):
	"""
	Takes an input filename and returns a string with the file extension removed.

	"""
	prevPos = 0
	currentPos = 0
	while currentPos != -1:
		prevPos = currentPos
		currentPos = inStr.find(".", prevPos+1)
	return inStr[0:prevPos]


if __name__ == "__main__":
	if len(sys.argv) == 2:
		main(sys.argv[1])
	else:
		print("Check usage: python cs_to_stars.py /path/to/your/cryosparc/file.cs")
		exit()

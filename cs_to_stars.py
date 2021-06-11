"""

Tristan Bell
Chao Lab
Massachusetts General Hospital

This script takes a cryosparc .cs file for particles or particle passthrough and exports the particle coordinates to a
group of star files that mimic a relion autopicking job.

"""

import sys
import numpy as np
import os
from os import listdir
from os.path import isfile, join


def main(inCs, inMcgs):
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

	# Read in the cs file as a np array
	print("Reading the cryosparc file...")
	f = np.load(inCs)

	# Parse the micrograph names from the ctf.star file
	print("Parsing micrograph names from star file...")
	mcg_parsed_names = parse_star(inMcgs)

	# Get the starting index for the particle location info
	print("Matching micrograph names between cs and star files...")
	start_index = infer_index(f)

	# Match star to cs entries and save the per-particle micrograph names
	mcg_names = []
	counter = 1
	previous_soln = "this_is_a_placeholder_name_hopefully_noone_ever_names_their_files_this_starfish_gorilla_massachusetts_taco_tarnish.taco"
	for i in range(0, len(f)):
		found_match_flag = False
		if counter % 100000 == 0:
			print(str(clean_large_numbers(counter))+" / "+str(clean_large_numbers(len(f))))
		# Try the previous solution before brute-forcing
		if no_ext(previous_soln) in str(f[i][start_index]):
			mcg_names.append(previous_soln)
			found_match_flag = True
		# Else brute-force it
		else:
			for mcg_name in mcg_parsed_names:
				# Save operations by doing calculations once
				target = str(f[i][start_index])
				if no_ext(mcg_name) in target:
					mcg_names.append(mcg_name)
					previous_soln = mcg_name
					found_match_flag = True
					break
		# Exit if the star and cs files don't have matching micrograph names
		if found_match_flag == False:
			print("\nERROR: There was no match found for the following micrograph: "+target)
			print("This can happen if your cryosparc and relion jobs didn't use the same input raw micrographs.")
			print("Perhaps the .cs file and the .star file don't correspond to the same data?")
			print("Exiting...")
			exit()
		counter += 1
	
	# Load the particles into a dictionary
	print("Transforming particle coordinates...")
	coord_dict = {}
	counter = 1
	for i in range(0, len(f)):
		# Pull micrograph name and add to dictionary
		mcg_name = mcg_names[i]
		if counter % 100000 == 0:
			print(str(clean_large_numbers(counter))+" / "+str(clean_large_numbers(len(f))))
		counter += 1
		if mcg_name not in coord_dict:
			coord_dict[mcg_name] = {"x":[], "y":[]}

		# Calculate transformed x and y coords
		h = f[i][start_index+1][0]
		l = f[i][start_index+1][1]
		x_frac = f[i][start_index+2]
		y_frac = f[i][start_index+3]
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
			print(clean_large_numbers(str(counter))+" / "+str(clean_large_numbers(len(coord_dict.keys()))))
		counter += 1
		g = open("Raw_data/"+no_ext(item)+"_autopick.star", "w", newline="")
		#g.write("\n# version 30001\n\ndata_particle\n\nloop_ \n_rlnCoordinateX #1 \n_rlnCoordinateY #2 \n_rlnAutopickFigureOfMerit #3 \n_rlnClassNumber #4 \n_rlnAnglePsi #5 \n")
		g.write("\n# version 30001\n\ndata_\n\nloop_ \n_rlnCoordinateX #1 \n_rlnCoordinateY #2 \n_rlnAutopickFigureOfMerit #3 \n_rlnClassNumber #4 \n_rlnAnglePsi #5 \n")
		for i in range(0, len(coord_dict[item]["x"])):
			g.write(line_writer(coord_dict[item]["x"][i], coord_dict[item]["y"][i]))
		g.write("\n")
		g.close()

	# Exit
	print("Done.")


def parse_star(inMcgs):
	# Open file
	f = open(inMcgs, "r")
	lines = f.readlines()

	# Find the entries starting as "MotionCorr"
	mcgs = []
	for i in range(0, len(lines)):
		if lines[i][:10] == "MotionCorr":
			# Split the line by " "
			temp_items = lines[i].split()
			mcgs.append(last_slash(temp_items[0]))
	return mcgs


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


def last_slash(inStr):
	prevPos = 0
	currentPos = 0
	while currentPos != -1:
		prevPos = currentPos
		currentPos = inStr.find("/", prevPos+1)
	return inStr[prevPos+1:]


def clean_large_numbers(inInt):
	"""
	Takes an integer and re-formats to string with human-readable comma-spaced numbers.

	"""
	inStr = str(inInt)
	outStr = ""
	
	if len(inStr) > 3:
		for i in range(1, len(inStr)+1):
			outStr = inStr[-i] + outStr
			if i % 3 == 0:
				outStr = "," + outStr
	else:
		outStr = inStr

	if outStr[0] == ",":
		return outStr[1:]
	else:
		return outStr


if __name__ == "__main__":
	if len(sys.argv) == 3:
		main(sys.argv[1], sys.argv[2])
	else:
		print("Check usage: python cs_to_stars.py /path/to/your/cryosparc/file.cs /path/to/your/micrographs_ctf.star")
		exit()

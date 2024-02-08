"""

Tristan Bell

This script takes a cryosparc .cs file for particles or particle passthrough and exports the particle coordinates to a
group of star files that mimic a relion autopicking job.

Revised to accommodate data format changes in RELION 4/5 and CryoSPARC 4.


"""

import sys
import numpy as np
import os
from os import listdir
from os.path import isfile, join
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--cs", help="cryosparc particle pick file")
parser.add_argument("--star", help="relion migrographs_ctf.star file")
parser.add_argument("--no_nudge", help="turn off dynamic selection in cs files", action="store_true")
parser.add_argument("--flipx", help="invert coordinates in x-axis", action="store_true")
parser.add_argument("--flipy", help="invert coordinates in y-axis", action="store_true")
parser.add_argument("--swapxy", help="swap x- and y-coordinates", action="store_true")
args = parser.parse_args()


def main(args):#inCs, inMcgs):#, args):
	# Parse arguments -> simple variables for back compatibility
	inCs = args.cs 
	inMcgs = args.star
	force_no_nudge = args.no_nudge

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
	start_index, x_index = infer_index(f)

	# Match star to cs entries and save the per-particle micrograph names
	mcg_names = []
	counter = 1
	previous_soln = "this_is_a_placeholder_name_hopefully_noone_ever_names_their_files_this_starfish_gorilla_massachusetts_taco_tarnish.mrc"
	for i in range(0, len(f)):
		#if counter < 20:
		found_match_flag = False
		if counter % 10000 == 0:
			print(str(clean_large_numbers(counter))+" / "+str(clean_large_numbers(len(f))))
		# Try the previous solution before brute-forcing
		if no_ext(previous_soln) in no_dot(str(f[i][start_index])):
			mcg_names.append(previous_soln)
			found_match_flag = True
		# Else brute-force it
		else:
			for mcg_name in mcg_parsed_names:
				# Save operations by doing calculations once
				target = no_dot(str(f[i][start_index]))
				if no_ext(mcg_name) in target:
					mcg_names.append(mcg_name)
					previous_soln = mcg_name
					found_match_flag = True
					break
		# Exit if the star and cs files don't have matching micrograph names
		if found_match_flag == False:
			try:
				print("\nERROR: There was no match found for the following micrograph: "+target)
			except:
				print("\nERROR: There was no match found for at least one micrograph specified in the cs file.")
			print("This can happen if your cryosparc and relion jobs didn't use the same input raw micrographs.")
			print("Perhaps the .cs file and the .star file don't correspond to the same data?")
			try:
				print("DEBUG: start_index = "+str(start_index))
			except:
				pass
			print("Exiting...")
			exit()
		counter += 1

	# Define indeces for x_frac and y_frac by checking for overly consistent coordinates
	# Check if x_index and x_index+1 are floats
	if (f[0][x_index].dtype == "float32") and (f[0][x_index] <= 1.0) and (f[0][x_index+1].dtype == "float32") and (f[0][x_index+1] <= 1.0):
		need_nudge_flag = True
	else:
		need_nudge_flag = False

	# If all three are coordinate floats of 1 or less, iterate over the particles and look repetition in positions +2, +3, and +4
	problem_indeces = []
	if start_index == x_index:
		if need_nudge_flag == True:
			for i in [2, 3, 4]:
				nonunique_items = []
				for j in range(0, len(f)):
					nonunique_items.append(f[j][start_index+i])
				if len(set(nonunique_items)) < len(f) * 0.01:
					problem_indeces.append(i)
			if force_no_nudge == True:
				print("\nWARNING: --no_nudge flag was invoked by user.\n\tThe script is forced to guess based on where they usually are.\n\tThe output coordinates may not be correct!\n\tPlease verify that your output pick positions match what you expect.\n\tProceeding...\n")
				adj_start_index = start_index
			elif (2 in problem_indeces) and (len(problem_indeces) == 1):
				print("\nNOTICE: Particle coordinates may have been stored in an unexpected place in the cryosparc file provided.\n\tWe have made our best guess as to where in the file your pick coordinates are.\n\tWe're quite confident they *should* be correct, but can't guarantee it.\n\tThis probably isn't a problem, but please carefully verify that your output pick positions match what you expect.\n\tTo disable this dynamic selection and fall back to where the coordinates *usually* are, rerun\n\t\t the program with the --no_nudge flag.\n\tIf this continues to be an issue, try a different cryosparc output file.\n\tProceeding...\n")
				adj_start_index = start_index + 1 # nudged start index to accomodate the extra cryosparc data
			else:
				print("\nWARNING: We could not unambiguously identify the particle coordinates in this cryosparc file.\n\tThe script is forced to guess based on where they usually are.\n\tThe output coordinates may not be correct!\n\tPlease verify that your output pick positions match what you expect.\n\tProceeding...\n")
				adj_start_index = start_index
		else:
			adj_start_index = start_index
	else:
		adj_start_index = x_index - 2
	
	# Load the particles into a dictionary
	print("Transforming particle coordinates...")
	coord_dict = {}
	counter = 1
	for i in range(0, len(f)):
		# Pull micrograph name and add to dictionary
		mcg_name = mcg_names[i]
		if counter % 50000 == 0:
			print(str(clean_large_numbers(counter))+" / "+str(clean_large_numbers(len(f))))
		counter += 1
		if mcg_name not in coord_dict:
			coord_dict[mcg_name] = {"x":[], "y":[]}

		# Calculate transformed x and y coords
		h = f[i][start_index+1][0]
		l = f[i][start_index+1][1]
		if args.flipx == False:
			x_frac = f[i][adj_start_index+2]
		else:
			x_frac = 1 - f[i][adj_start_index+2]
		if args.flipy == False:
			y_frac = f[i][adj_start_index+3]
		else:
			y_frac = 1 - f[i][adj_start_index+3]

		if args.swapxy == False: #swapxy implementation
			x_coord = round((l*x_frac), 0)
			y_coord = round(h-(h*y_frac), 0)
		else:
			x_coord = round(h-(h*y_frac), 0)
			y_coord = round((l*x_frac), 0)

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

	# Find the position of "_rlnMicrographName" in the star file loop definition
	for i in range(0, len(lines)):
		if "loop_" in lines[i]:
			loop_start = i
			for j in range(i+1, len(lines)):
				if "_rlnMicrographName " in lines[j]:
					parse_pos = int(lines[j][lines[j].find("#")+1:])-1
					break

	# Isolate the micrograph name from each line of the star file
	mcgs = []
	for i in range(loop_start+1, len(lines)):
		if lines[i][0] != "_":
			mcgs.append(last_slash(lines[i].split(" ")[parse_pos]))

	# Final check for zero-length strings
	clean_mcgs = []
	for i in range(0, len(mcgs)):
		if len(mcgs[i]) > 0:
			clean_mcgs.append(mcgs[i])
	
	return clean_mcgs


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


def assess_variability(np_array, start_index):
	# Check the variability in floats at positions i+2, i+3, i+4
	# If i+2 and i+3 are variable but i+4 is not, return False
	# Else return True
	plus2 = []
	plus3 = []
	plus4 = []

	for i in range(len(np_array)):
		plus2.append(np_array[i][start_index+2])
		plus3.append(np_array[i][start_index+3])
		plus4.append(np_array[i][start_index+4])

	plus2setlen = len(set(plus2))
	plus3setlen = len(set(plus3))
	plus4setlen = len(set(plus4))

	if plus2setlen > plus4setlen:
		return False
	else:
		return True
	

def infer_index(np_array):
	# Default defined pattern is binary string, list of two ints >1000, float <= 1, float <=1
	# If there's a subsequent float 0-1, then scan the array to see which sites are most variable
	for i in range(0, len(np_array[1])):
		try:
			np_array[1][i].decode("utf-8")
			try:
				a = len(np_array[1][i+1])
				if (np_array[1][i+2] >= 0) and (np_array[1][i+2] <= 1):
					if (np_array[1][i+3] >= 0) and (np_array[1][i+3] <= 1):
						startInd = i
						if (np_array[1][i+4] >= 0) and (np_array[1][i+4] <= 1):
							if assess_variability(np_array, i) == False:
								return i, i+2
							else:
								return i, i+3
						break
			except:
				pass
		except:
			pass
	try:
		return startInd, startInd+2
	except: # backup pattern is binary string, list of two ints, float > 1, float <1, float <1
		for i in range(0, len(np_array[1])):
			try:
				np_array[1][i].decode("utf-8")
				try:
					a = len(np_array[1][i+1])
					if np_array[1][i+2] >= 1:
						if (np_array[1][i+3] <= 1) and (np_array[1][i+3] >= 0):
							if (np_array[1][i+4] <= 1) and (np_array[1][i+4] >= 0):
								if assess_variability(np_array, i) == False:
									return i, i+3
								else:
									return i, i+4
				except:
					pass
			except:
				pass


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


def no_dot(inStr):
	"""
	Relion converts "." in cryoparc names to "_" - this function takes a script and performs this
	conversion prior to name-matching.

	"""
	return inStr.replace(".", "_")


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
	"""
	Returns the component of a string past the last forward slash character.
	"""
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
	if (args.cs == None) or (args.star == None):
		print("Check usage: python cs_to_stars.py --cs /path/to/your/cryosparc/file.cs --star /path/to/your/micrographs_ctf.star\nUse python cs_to_stars.py --help for all options.")
	else:
		main(args)

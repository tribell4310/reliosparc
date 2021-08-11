"""

converge.py

Script for tracking convergence between classes in relion in a continue run.

"""


import sys
from os.path import isfile, join
import time
import matplotlib as mpl
mpl.use("Agg")
mpl.rc('figure', max_open_warning = 0)
import matplotlib.pyplot as plt



def main(prevIterations, totalIterations):
	# Set a counter of iterations and an empty dictionary
	iteration = 1
	newDict = {}

	# Create a logfile
	g = open("converge.log", "w")
	g.write(" ")
	g.close()

	# Set up a loop that loads each data file as it appears
	while iteration <= prevIterations:
		current_star = "run_it"+threepad(str(iteration), 3)+"_data.star"
		if isfile("run_it"+threepad(str(iteration+1), 3)+"_data.star") == False:
			time.sleep(20)

		if isfile(current_star):
			# Pause for 15s to allow relion to actually write the star file.
			if isfile("run_it"+threepad(str(iteration+1), 3)+"_data.star") == False:
				time.sleep(15)

			# Read the current star file and populate the dictionary
			newDict = parse_star(current_star, newDict)
			
			# Calculate the number of particles that switched classes in the last iteration
			change_counter = 0
			total_counter = 0
			if iteration > 1:
				for item in newDict:
					total_counter += 1
					if newDict[item][iteration-1] != newDict[item][iteration-2]:
						change_counter += 1
				g = open("converge.log", "a")
				g.write("Iteration: "+str(iteration)+"\n")
				g.write(str(change_counter)+" of "+str(total_counter)+" particles were re-assigned to new classes.\nThis is "+str(round(100*float(change_counter)/float(total_counter), 1))+" percent of the input pool.\n\n")
				g.close()
				# Count the number of particles in each class
				counts_dict = class_counter(counts_dict, newDict)
			else:
				for item in newDict:
					total_counter += 1
				g = open("converge.log", "a")
				g.write("Iteration: "+str(iteration)+"\n")
				g.write("Loaded data for "+str(total_counter)+" particles.\n\n")
				g.close()
				# Count the number of particles in each class
				counts_dict = class_counter({}, newDict)
			
			# Create a plot of the class populations over time
			x = range(1, iteration+1)
			fig, ax = plt.subplots()
			ax.set(xlabel="Iteration", ylabel="Particles in Class")
			ax.grid()
			for item in counts_dict:
				y = counts_dict[item]
				ax.plot(x, y, label=item)
			ax.legend()
			fig.savefig("converge_"+threepad(str(iteration), 3)+".png")
			plt.clf()

			#Increment the iteration counter 
			iteration += 1

	# Load and process the continue iterations
	while iteration <= totalIterations:
		current_star = "run_ct"+str(prevIterations)+"_it"+threepad(str(iteration), 3)+"_data.star"
		if isfile("run_ct"+str(prevIterations)+"_it"+threepad(str(iteration+1), 3)+"_data.star") == False:
			time.sleep(20)

		if isfile(current_star):
			# Pause for 15s to allow relion to actually write the star file.
			if isfile("run_ct"+str(prevIterations)+"_it"+threepad(str(iteration+1), 3)+"_data.star") == False:
				time.sleep(15)

			# Read the current star file and populate the dictionary
			newDict = parse_star(current_star, newDict)
			
			# Calculate the number of particles that switched classes in the last iteration
			change_counter = 0
			total_counter = 0
			if iteration > 1:
				for item in newDict:
					total_counter += 1
					if newDict[item][iteration-1] != newDict[item][iteration-2]:
						change_counter += 1
				g = open("converge.log", "a")
				g.write("Iteration: "+str(iteration)+"\n")
				g.write(str(change_counter)+" of "+str(total_counter)+" particles were re-assigned to new classes.\nThis is "+str(round(100*float(change_counter)/float(total_counter), 1))+" percent of the input pool.\n\n")
				g.close()
				# Count the number of particles in each class
				counts_dict = class_counter(counts_dict, newDict)
			
			# Create a plot of the class populations over time
			x = range(1, iteration+1)
			fig, ax = plt.subplots()
			ax.set(xlabel="Iteration", ylabel="Particles in Class")
			ax.grid()
			for item in counts_dict:
				y = counts_dict[item]
				ax.plot(x, y, label=item)
			ax.legend()
			fig.savefig("converge_"+threepad(str(iteration), 3)+".png")
			plt.clf()

			#Increment the iteration counter 
			iteration += 1


def class_counter(countDict, inDict):
	# Is the countDict empty?  Is so, populate it with empty lists.
	if countDict == {}:
		newDict = {}
		for item in inDict:
			if inDict[item][-1] not in newDict.keys():
				newDict[inDict[item][-1]] = []
	else:
		newDict = countDict.copy()

	# Add a zero count to each item in newDict
	for item in newDict:
		newDict[item].append(0)

	# Count the current class populations
	for item in inDict:
		newDict[inDict[item][-1]][-1] += 1

	return newDict

def parse_star(inStar, prevDict):
	# Open file
	f = open(inStar, "r")
	lines = f.readlines()

	# Find the position of "_rlnClassNumber" and "_rlnImageName" in the star file loop definition
	image_found = False
	class_found = False
	for j in range(0, len(lines)):
		if "_rlnImageName" in lines[j]:
			if image_found == False:
				image_parse_pos = int(lines[j][lines[j].find("#")+1:])-1
				image_index = j
				image_found = True
		if "_rlnClassNumber" in lines[j]:
			if class_found == False:
				class_parse_pos = int(lines[j][lines[j].find("#")+1:])-1
				class_index = j
				class_found = True
		if (image_found == True) and (class_found == True):
			break

	# Is the dictionary empty?  If so, populate it with the micrograh names.
	if prevDict == {}:
		for i in range(max(image_index, class_index)+1, len(lines)):
			if (lines[i][0] != "_") and len(lines[i]) > 20:
				prevDict[lines[i].split()[image_parse_pos]] = []

	# Isolate the features from each line of the star file
	for i in range(max(image_index, class_index)+1, len(lines)):
			if (lines[i][0] != "_") and len(lines[i]) > 20:
				this_class = lines[i].split()[class_parse_pos]
				prevDict[lines[i].split()[image_parse_pos]].append(this_class)

	return prevDict


def threepad(inStr, final_len):
	while len(inStr) < final_len:
		inStr = "0"+inStr
	return inStr


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
	if len(sys.argv) == 3:
		print("Live output will be logged in: converge.log")
		main(int(sys.argv[1]), int(sys.argv[2]))
	else:
		print("Check usage: python converge_continue.py numberOfInitialIterations numberOfTotalIterations")
		exit()

# Python code that is meant as my program, in the python course I am currently taking


# This program takes a list of Protein names, as well as their abundancies.
# It then compares protein levels, between different samples.

# The program should have an option to choose between analysing simulation data or experimental data
# This could for example tell the program if normalization is required or not.

file = '/Users/johanbj/Documents/Courses/Python/protein_test.csv'

import csv
with open(file, newline='') as csvfile:
	spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')
	for row in spamreader:
		print(row[0])
		
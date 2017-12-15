# Python code that is meant as my program, in the python course I am currently taking

# This program invstigates the different pathways given in the GEM "ecYeast_v.1_4" (GECKO) model.
# It the gatters and compares som data on the different enzymes which are annotated to each pathway.
#
# Data to compare, which is given in the model, include; molecular weights, length of enzymes and amino
# acid usage. 



############# READ THE TEXT FILE CONTAINING INFORMATION ON PATHWAYS THAT PROTEINS PARTICIPATE IN #############

read_this_file = '../Protein_Pathways.txt'

protein_pathways = open(read_this_file,'r')

for line in protein_pathways:     # Go through all lines in the text-file
	line = line.rstrip()          # Get rid of the new-line character.

	#print(repr(line))
	#print(repr(line.rstrip()))
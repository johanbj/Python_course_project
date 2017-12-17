# Python code that is meant as my program, in the python course I am currently taking

# This program invstigates the different pathways given in the GEM "ecYeast_v.1_4" (GECKO) model.
# It the gatters and compares som data on the different enzymes which are annotated to each pathway.
#
# Data to compare, which is given in the model, include; molecular weights, length of enzymes and amino
# acid usage. 



############# READ THE TEXT FILE CONTAINING INFORMATION ON PATHWAYS THAT PROTEINS PARTICIPATE IN #############

read_this_file = '../Protein_Pathways.txt'

protein_pathways = open(read_this_file,'r')
collection_of_pathways = []

for line in protein_pathways:     # Go through all lines in the text-file
	line = line.rstrip()          # Get rid of the new-line character.
	list_from_line = line.split(' sce0')

#############################################

# I NEED TO ADD THE 'SCE' TO THE ELEMENT AGAIN

#############################################
	
	collection_of_pathways.extend(list_from_line)
print(collection_of_pathways)
print(len(collection_of_pathways))

unique_pathways = []
for pathway_code in collection_of_pathways:
    if pathway_code not in unique_pathways:
        unique_pathways.append(pathway_code)
print(unique_pathways)
print(len(unique_pathways))
	#for pathway in list_from_line:
	#	if pathway
	#print(repr(line))
	#print(repr(line.rstrip()))



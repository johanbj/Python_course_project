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
	list_from_line = line.split(' sce0')     # Split and categorize the line, based on a unique identifier

	if listlen>1:                       # If the list now has more than 1 element, all elements except
		for j in range(1,(listlen)):    # the first must have 'sce0' added to it again.
			list_from_line[j] = 'sce0' + list_from_line[j]
	
	collection_of_pathways.extend(list_from_line)    # Collect all split pathway names

unique_pathways = []
for pathway_name in collection_of_pathways:
    if pathway_name not in unique_pathways:     # Filter so that all names will be unique in the list
        unique_pathways.append(pathway_name)

protein_pathways.close()     # Close the file
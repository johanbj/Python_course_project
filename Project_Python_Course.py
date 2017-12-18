# Python code that is meant as my program, in the python course I am currently taking

# This program invstigates the different pathways given in the GEM "ecYeast_v.1_4" (GECKO) model.
# It the gatters and compares som data on the different enzymes which are annotated to each pathway.
#
# Data to compare, which is given in the model, include; molecular weights, length of enzymes and amino
# acid usage. 



############  READ THE TEXT FILE CONTAINING INFORMATION ON PATHWAYS THAT PROTEINS PARTICIPATE IN  ############

read_the_pathways = '../Protein_Pathways.txt'

protein_pathways = open(read_the_pathways,'r')
collection_of_pathways = []
collection_of_lines = []

for line in protein_pathways:          # Go through all lines in the text-file
	line = line.rstrip()               # Get rid of the new-line character.
	collection_of_lines.append(line)   # Save the original lines for later purpose
	list_from_line = line.split(' sce0')     # Split and categorize the line, based on a unique identifier
	listlen = len(list_from_line)

	if listlen>1:                       # If the list now has more than 1 element, all elements except
		for j in range(1,(listlen)):    # the first must have 'sce0' added to it again.
			list_from_line[j] = 'sce0' + list_from_line[j]
	
	collection_of_pathways.extend(list_from_line)    # Collect all split pathway names

unique_pathways = []
for pathway_name in collection_of_pathways:
    if pathway_name not in unique_pathways:     # Filter so that all names will be unique in the list
        unique_pathways.append(pathway_name)

protein_pathways.close()     # Close the file

###############################################################################################################

# Next task is to go through all the unique pathway names & note which protein that
# are associated to that name.


##########  MATCH AND ORDER PROTEIN-INDEX CORRESPONDING TO EACH PATHWAY  ##########

Pathways_and_enzymes = {}

for p_name in unique_pathways:
	i = 1
	hits = []
	for enzyme in collection_of_lines:
		if p_name in enzyme:
			hits.append(i)

		i += 1
	Pathways_and_enzymes[p_name] = hits

###################################################################################

# The next part of the analysis is to compute interesting data for each pathway
# One example may be the average enzyme molecular weight, for a pathway
# To begin with, code will be tested with only one pathway. The pathway was
# Choosen by random and is 'sce04011  MAPK signaling pathway - yeast'.



##########  CALCULATE THE AVERAGE ENZYME MOLECULAR WEIGHT OF THE PATHWAY  ##########

####  EXTRACT MOLECULAR WEIGHTS OF ENZYMES  ####
read_the_MWs = '../Protein_MWs.txt'         
protein_MWs = open(read_the_MWs,'r')
collected_weights = []

for weight in protein_MWs:
	weight = weight.rstrip()
	weight = float(weight)
	collected_weights.append(weight)

protein_MWs.close()
########################################

########  EXTRACT SEQUENCES OF ENZYMES  ########
read_the_seqs = '../Protein_Sequences.txt'         
protein_seqs = open(read_the_seqs,'r')
collected_seqs = []

for seq in protein_seqs:
	seq = seq.rstrip()
	collected_seqs.append(seq)

protein_seqs.close()
################################################

###################################################################################

# Now when the information about all enzymes has been provided, the script
# will calculate the data for each pathway and its connected enzymes


################  CALCULATING DATA  ################

import numpy

def averages_of_enzymes(list_of_enzymes):      # Defining a function that will do ALL calculations
	total_weight = 0																					# I
	total_length = 0																					# N
	total_amount = 0																					# I
	list_of_amino_acids = [0] * 20																		# T
	amino_acids = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']		# !

	for ez in list_of_enzymes:	# Goes through all enzymes connected to the pathway investigated
		ez = int(ez)									
		total_weight += collected_weights[(ez-1)]       # Add weight to the total weight
		total_length += len(collected_seqs[(ez-1)])     # Add length to the total length

		for letter in amino_acids:												#
			abundancy_of_amino_acid = collected_seqs[(ez-1)].count(letter)		#
			INDEX = amino_acids.index(letter)									#
			list_of_amino_acids[INDEX] += abundancy_of_amino_acid				#

	list_of_amino_acids = numpy.array(list_of_amino_acids)
	list_of_amino_acids = (list_of_amino_acids/total_length)*100
	average_ez_weight = total_weight/len(list_of_enzymes)
	average_ez_length = total_length/len(list_of_enzymes)


	averages_of_enz = [len(list_of_enzymes),average_ez_weight,average_ez_length]
	averages_of_enz.extend(list_of_amino_acids)
	return(averages_of_enz)


used_enzymes = Pathways_and_enzymes['sce04011  MAPK signaling pathway - yeast']

averages_of_enz = averages_of_enzymes(used_enzymes)




Pathways_and_enzymes_analyzed = {}
# The dictionary contain data on each pathway, in the following order; number of enzymes in pathway,
# average enzyme weight, average enzyme length percentage of different amino acids needed in the order
# A, R, N, D, C, E, Q, G, H, I, L, K, M, F, P, S, T, W, Y, V




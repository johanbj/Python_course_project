# Python code that is meant as my program, in the python course I am currently taking

# This program investigates the different pathways given in the GEM "ecYeast_v.1_4" (GECKO) model.
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

unique_pathways = unique_pathways[:-1]          # Remove the last element which is a blank character


###############################################################################################################

# Next task is to go through all the unique pathway names & note which protein that
# are associated to that name.


##########  MATCH AND ORDER PROTEIN-INDEX CORRESPONDING TO EACH PATHWAY  ##########


Pathways_and_enzymes = {}
cutoff_enzymes_per_pathway = 1
save_these_pathways = []
pathway_index = 0

for p_name in unique_pathways:               # Go through each unique pathway name
	i = 1
	hits = []
	for enzyme in collection_of_lines:		 # Check for each enzyme if they are present in the pathway
		if p_name in enzyme:
			hits.append(i)					 # If yes, then store the enzyme index (first enzyme has index 1)

		i += 1

	if len(hits) >= cutoff_enzymes_per_pathway:
		Pathways_and_enzymes[p_name] = hits
		save_these_pathways.append(pathway_index)
	
	pathway_index += 1

updated_unique_pathways = []

for pathway_to_save in save_these_pathways:
	updated_unique_pathways.append(unique_pathways[pathway_to_save])
	
unique_pathways = updated_unique_pathways


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

########  EXTRACT GENENAMES OF ENZYMES  ########
read_the_genenames = '../Protein_GeneNames.txt'         
protein_genenames = open(read_the_genenames,'r')
collected_genenames = []

for genename in protein_genenames:
	genename = genename.rstrip()
	collected_genenames.append(genename)

protein_genenames.close()
################################################


###################################################################################


# Now when the information about all enzymes has been provided, the script
# will calculate the data for each pathway and its connected enzymes


################  CALCULATING DATA  ################


import numpy


######  FUNCTION THAT COLLECTS ALL INTERESTING DATA


def averages_of_enzymes(list_of_enzymes):      # Defining a function that will do ALL calculations
	total_weight = 0																					# I
	total_length = 0																					# N
	total_amount = 0																					# I
	list_of_amino_acids = [0] * 20																		# T
	amino_acids = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']		# !

	for ez in list_of_enzymes:	# Goes through all enzymes connected to the pathway investigated
		ez = int(ez)									
		total_weight += collected_weights[(ez-1)]       # Add weight to the total weight (ez-1 to get the correct index)
		total_length += len(collected_seqs[(ez-1)])     # Add length to the total length (ez-1 to get the correct index)

		for letter in amino_acids:												# For each different amino aicd
			abundancy_of_amino_acid = collected_seqs[(ez-1)].count(letter)		# Count its presence in the indexed enzyme
			INDEX = amino_acids.index(letter)									# Find the index of the investigated amino acid
			list_of_amino_acids[INDEX] += abundancy_of_amino_acid				# Save the acquired value to that index

	list_of_amino_acids = numpy.array(list_of_amino_acids)             # 
	list_of_amino_acids = (list_of_amino_acids/total_length)*100       # Normalize amino acid abundacies to percentage
	average_ez_weight = total_weight/len(list_of_enzymes)              # Get the average weight of an enzyme in the pathway
	average_ez_length = total_length/len(list_of_enzymes)              # Get the average length of an enzyme in the pathway


	averages_of_enz = [len(list_of_enzymes),average_ez_weight,average_ez_length]    # Collect all data
	averages_of_enz.extend(list_of_amino_acids)										# Collect all data
	return(averages_of_enz)


####################################################


# As the function for calculating the data, that I am interested in, is complete I go through all the
# pathways and compute it for them. I store the information in the dictionary initiated below.
# The dictionary will contain data on each pathway, in the following order; number of enzymes in pathway,
# average enzyme weight, average enzyme length percentage of different amino acids needed in the order
# A, R, N, D, C, E, Q, G, H, I, L, K, M, F, P, S, T, W, Y, V


Pathways_and_enzymes_analyzed = {}


for thepathway in unique_pathways:
	used_enzymes = Pathways_and_enzymes[thepathway]
	averages_of_enz = averages_of_enzymes(used_enzymes)
	Pathways_and_enzymes_analyzed[thepathway]=averages_of_enz


#########################################################################################
#########################################################################################


# After having computed all the relevant data, I now want to do some
# basic plotting in order to visualize my results.

import matplotlib.pyplot as plt


###################  PLOT AMOUNT OF PROTEINS PRESENT FOR ALL PATHWAYS  ###################


number_of_proteins = []
weights_of_proteins = []
lengths_of_proteins = []
number_times_lengths_of_proteins = []

for thepathway in unique_pathways:
	data_from_pathway = Pathways_and_enzymes_analyzed[thepathway]
	number_of_proteins.append(data_from_pathway[0])
	weights_of_proteins.append(data_from_pathway[1])
	lengths_of_proteins.append(data_from_pathway[2])
	number_times_lengths = data_from_pathway[0]*data_from_pathway[2]
	number_times_lengths_of_proteins.append(number_times_lengths)

#number_of_proteins = number_of_proteins[:-1]
#weights_of_proteins = weights_of_proteins[:-1]
#lengths_of_proteins = lengths_of_proteins[:-1]
#number_times_lengths_of_proteins = number_times_lengths_of_proteins[:-1]


data_tags = ['Number_of_proteins','Weights_of_proteins','Lengths_of_proteins','Number_times_lengths_of_proteins']
data_to_plot = {'Number_of_proteins': number_of_proteins,
				'Weights_of_proteins': weights_of_proteins,
				'Lengths_of_proteins': lengths_of_proteins,
				'Number_times_lengths_of_proteins': number_times_lengths_of_proteins}

for tag in data_tags:
	pathway_pos = numpy.arange(1,len(unique_pathways)+1)
	plt.bar(pathway_pos, data_to_plot[tag], align='center', alpha=0.5)
	plt.ylabel(tag)
	title_line = tag + ' in pathways'
	plt.title(title_line)
	plt.grid(axis='y')

	orig_high_abound = max(data_to_plot[tag])
	high_abound = round(orig_high_abound)

	if high_abound < orig_high_abound:
		high_abound += 1

	plt.axis([0, (len(unique_pathways)+1), 0, high_abound])
	save_name = tag + '_in_pathways.png'
	plt.savefig(save_name, bbox_inches='tight')
	plt.close()


###########################  CHECK THE 10 MOST EXPENSIVE PATHWAYS  #############################


sort_index = numpy.argsort(number_times_lengths_of_proteins)
printthis = sort_index[:-11:-1]

for thisnumber in printthis:
	print(unique_pathways[thisnumber])


#####################  CHECK AMINO ACID ABUNDANCY IN ALL PROTEINS AT INCREASING SIZE  #####################


amino_acids = ('A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')

lengths_of_proteins_check = []

for seq_length in collected_seqs:
	lengths_of_proteins_check.append(len(seq_length))

sort_lengths_of_proteins = numpy.argsort(lengths_of_proteins_check)
sorted_lengths_of_proteins = sorted(lengths_of_proteins_check)


########  Plotting the different protein lengths


prot_number = len(lengths_of_proteins_check) + 5
pathway_pos = numpy.arange(5,prot_number)
plt.plot(pathway_pos, sorted_lengths_of_proteins)
plt.ylabel('Protein lengths')
title_string = 'Protein distribution'
plt.title(title_string)
plt.grid(axis='y')

orig_high_abound = max(lengths_of_proteins_check)
high_abound = round(orig_high_abound)

if high_abound < orig_high_abound:
	high_abound += 1

prot_bound = prot_number + 10
plt.axis([0, prot_bound, 0, high_abound])
name_of_fig = title_string + '.png'
plt.savefig(name_of_fig, bbox_inches='tight')
plt.close()


##################################################


amino_acids_in_all_proteins = {}

for letter in amino_acids:
	list_of_amino_acids_protein_length = []
	for sorted_length in sort_lengths_of_proteins:												
		abundancy_of_amino_acid_protein_length = collected_seqs[sorted_length].count(letter)		# Count its presence in the indexed enzyme
		abundancy_of_amino_acid_protein_length = (abundancy_of_amino_acid_protein_length/len(collected_seqs[sorted_length]))*100								
		list_of_amino_acids_protein_length.append(abundancy_of_amino_acid_protein_length)			# Save the acquired value to that index
	amino_acids_in_all_proteins[letter] = list_of_amino_acids_protein_length


##########  PLOT THE RESULTS FOR EACH AMINO ACID & PLOT AS A BOXPLOT

boxplot_dataset_aa_in_proteins = []

# Regular plot


from scipy import stats

for aa in amino_acids:
	slope, intercept, r_value, p_value, std_err = stats.linregress(pathway_pos,amino_acids_in_all_proteins[aa])
	prot_number = len(lengths_of_proteins_check) + 5
	pathway_pos = numpy.arange(5,prot_number)
	amino_acid_abundancy = amino_acids_in_all_proteins[aa]
	fit = numpy.polyfit(pathway_pos,amino_acids_in_all_proteins[aa],1)
	fit_fn = numpy.poly1d(fit)
	plt.scatter(pathway_pos, amino_acids_in_all_proteins[aa],s=1,c='blue',marker='.')
	plt.plot(pathway_pos,fit_fn(pathway_pos),c='black')
	text_in_plot = 'r-squared: ' + str(r_value**2)
	plt.text(10, 10, text_in_plot, fontsize=10)
	#plt.bar(pathway_pos, amino_acids_in_all_proteins[aa], align='center', alpha=0.5)
	plt.ylabel('Average Amino Acid aundancy in %')
	title_string = 'Amino_acid_' + aa + '_in proteins'
	plt.title(title_string)
	plt.grid(axis='y')

	orig_high_abound = max(amino_acid_abundancy)
	high_abound = round(orig_high_abound)

	if high_abound < orig_high_abound:
		high_abound += 1

	prot_bound = prot_number + 10
	plt.axis([0, prot_bound, 0, high_abound])
	name_of_fig = title_string + '.png'
	plt.savefig(name_of_fig, bbox_inches='tight')
	plt.close()

	print_statement = 'Done plotting amino acid ' + aa + ' for proteins'
	print(print_statement)

	boxplot_dataset_aa_in_proteins.append(amino_acids_in_all_proteins[aa])


# Boxplot

flierprops = dict(marker='.', markerfacecolor='black', markersize=1,linestyle='none')	
plt.boxplot(boxplot_dataset_aa_in_proteins,flierprops=flierprops,widths=0.75)
plt.xlim([0,21])
amino_acids_pos = numpy.arange(1,len(amino_acids)+1)
plt.xticks(amino_acids_pos, amino_acids)
plt.grid(axis='y')
plt.ylabel('Abundance in protein (%)')
plt.title('Box plot representation of amino acid abundancies in proteins')
plt.savefig('Box_plot_of_amino_acid_abundancies_proteins.png', bbox_inches='tight')
plt.close()


###########################################################################################################



# Then check the amino acids for each pathway

################  PLOT ABUNDANCE OF A GIVEN AMINO ACID FOR ALL PATHWAYS  #################
amino_acids = ('A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
results_amino_acid_for_boxplot = {}

for aa in amino_acids:
	results_amino_acid_for_boxplot[aa] = []


for aa in amino_acids:
	investigate_this_amino_acid = aa

	# investigate_this_amino_acid = 'W'
	# Input one of the above letters

	INDEX = amino_acids.index(investigate_this_amino_acid)      # This feature is added because it is needed if the
															    # user would like to input a specific amino acid

	investigated_amino_acid_abundancy = []

	for thepathway in unique_pathways:
		data_from_pathway = Pathways_and_enzymes_analyzed[thepathway]
		data_from_pathway = data_from_pathway[3:]

		for aa_index in range(0,20,1):
			old_values = results_amino_acid_for_boxplot[amino_acids[aa_index]]
			old_values.append(data_from_pathway[aa_index])
			results_amino_acid_for_boxplot[amino_acids[aa_index]] = old_values

		investigated_amino_acid_abundancy.append(data_from_pathway[INDEX])


	pathway_pos = numpy.arange(1,len(unique_pathways)+1)
	amino_acid_abundancy = investigated_amino_acid_abundancy
	plt.bar(pathway_pos, amino_acid_abundancy, align='center', alpha=0.5)
	plt.ylabel('Average Amino Acid aundancy in %')
	title_string = 'Amino_acid_' + investigate_this_amino_acid + '_in pathways'
	plt.title(title_string)
	plt.grid(axis='y')

	orig_high_abound = max(investigated_amino_acid_abundancy)
	high_abound = round(orig_high_abound)

	if high_abound < orig_high_abound:
		high_abound += 1

	plt.axis([0, (len(unique_pathways)+1), 0, high_abound])
	name_of_fig = title_string + '.png'
	plt.savefig(name_of_fig, bbox_inches='tight')
	plt.close()

	print_statement = 'Done plotting amino acid ' + aa + ' for pathways'
	print(print_statement)


################  CREATE A BOXPLOT FOR ABUNDANCIES OF ALL AMINO ACIDS IN THE PATHWAYS  #################

boxplot_dataset = []
for aa in amino_acids:
	boxplot_dataset.append(results_amino_acid_for_boxplot[aa])


flierprops = dict(marker='.', markerfacecolor='black', markersize=1,linestyle='none')	
plt.boxplot(boxplot_dataset,flierprops=flierprops,widths=0.75)
plt.xlim([0,21])
amino_acids_pos = numpy.arange(1,len(amino_acids)+1)
amino_acid_abundancy = boxplot_dataset
plt.xticks(amino_acids_pos, amino_acids)
plt.grid(axis='y')
plt.ylabel('Abundance in pathway (%)')
plt.title('Box plot representation of amino acid abundancies in pathways')
plt.savefig('Box_plot_of_amino_acid_abundancies_pathways.png', bbox_inches='tight')
plt.close()


##########################################################################################



########################################
#####   EXAMPLE ANALYSIS OF DATA   #####
########################################


if cutoff_enzymes_per_pathway == 5:

# After having a look at the plotted data I will further examine some features


# I would like to view the protein that has the highest abundancy for each
# amino acid. Since each amino acid, by looking at the boxplot, seem to have
# some high "outlier"

# Isoleucine, Leucine and Phenylalanine all have peak abundancies in pathway
# number 64, therefore I will look closer into this one.

# I would like to investigate the pathway that only has 3 % of Aspartate

##################
#######  High abundancy proteins

#######  Sort the protein sequences so that they are in order of ascending length
	sorted_sequences = numpy.argsort(lengths_of_proteins_check)
	sorted_collected_seqs = []
	sorted_collected_genenames = []

	for number_from_sorted_sequences in sorted_sequences:
		sorted_collected_seqs.append(collected_seqs[number_from_sorted_sequences])
		sorted_collected_genenames.append(collected_genenames[number_from_sorted_sequences])
#################################################################################

#####  Print the protein sequence of the enzyme that contains the highest amount of each amino acid #####
	amino_acids = ('A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
	for aa in amino_acids:
		sorted_abundancies_for_amino_acid = numpy.argsort(amino_acids_in_all_proteins[aa])
		print(' ')
		print('The protein that has the highest percentual amount of ',aa, ' is: ',sorted_collected_genenames[sorted_abundancies_for_amino_acid[-1]])
		print(' ')
		print('Which has the following amino acid sequence:')
		print(' ')
		print(sorted_collected_seqs[sorted_abundancies_for_amino_acid[-1]])
#########################################################################################################

##################


###############
#######  Pathway 64 (Isoleucine, Leucine and Phenylalanine)

# I investigate the seemingly interesting pathway 64

	check_this_pathway = 64

	print(' ')
	print('The pathway high in Isoleucine, Leucine and Phenylalanine is:')
	print(' ')
	print(unique_pathways[(check_this_pathway-1)])
	print(' ')
	print('This pathway concists of the proteins:')
	print(' ')

	enzymes_in_this_pathway = Pathways_and_enzymes[unique_pathways[(check_this_pathway-1)]]

	for enzyme_number in enzymes_in_this_pathway:
		print(collected_genenames[enzyme_number-1])
		print(' ')

###############


###############
#######  Aspartate
# Instead of searching in which pathway Aspartate has its lowest abundancy I choose instead to
# look at the plot generated by this script, which shows that I should investigate pathway 65

	check_this_pathway_aspartate = 65

	print(' ')
	print('The pathway low in Aspartate is:')
	print(' ')
	print(unique_pathways[(check_this_pathway_aspartate-1)])
	print(' ')
	print('This pathway concists of the proteins:')
	print(' ')

	enzymes_in_this_pathway = Pathways_and_enzymes[unique_pathways[(check_this_pathway_aspartate-1)]]

	for enzyme_number in enzymes_in_this_pathway:
		print(collected_genenames[enzyme_number-1])
		print(' ')

###############









import csv
import json 
import os
import subprocess
from Bio import SeqIO
import copy
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

from func3_quantification import inputFileHandle

def func4_outputResult(inputPeptide, inputCsvFile, tumourResultDirectories, normalResultDirectories, outputPath, withZero, dpi, selectedPeptide, selectedRegion):
	# STEP 1 : get the input peptide sequence
	TSA_list = getTSAlist(inputPeptide)

	# STEP 2 : input the sample result directories file and load the output dictionary for that specific sample
	tumour_result = extractResult(tumourResultDirectories, TSA_list)
	normal_result = extractResult(normalResultDirectories, TSA_list)
	total_result = combine_tumour_normal_results(tumour_result, normal_result)

	
	# STEP 3 : get the padding sequence for all region in the result and output as fasta files for each region
	get_padding(tumour_result, outputPath)
	
	# STEP 4.1 : create the overall csv
	df1 = output_overall_csv(tumour_result, total_result,TSA_list, outputPath, withZero)
	# STEP 4.2 : create the peptide maximum region csv file
	df2 = output_total_max_csv(tumour_result, total_result, TSA_list, outputPath)
	# STEP 4.3 : create the sample maximum region csv file
	df3 = output_sample_max_csv(total_result, outputPath)

	# STEP 5.0 : preprocess data for the plots
	try:
		inputCsvResult, inputCsvError = inputFileHandle(inputCsvFile)
		if (inputCsvError is not None):
			raise ValueError(f"input CSV file error : {inputCsvError}")
	except ValueError as err:
		print(f"{err}", file = sys.stderr)
		return err 
	except Exception:
		print(f"An error occurred, please check your inputs", file = sys.stderr)
		return "An error occurred, please check your inputs"
		 
	tumour = inputCsvResult[0]
	normal = inputCsvResult[1]
	
	# STEP 5.1 : create plot 1
	dpi = int(dpi)
	withZero = int(withZero)
	try:
		if (selectedRegion != None and selectedPeptide == None):
			raise ValueError("Argument error: missing selectedPeptide, please include your selected peptide for the selected regions")
		
		if (selectedPeptide != None):
			selectedRegion = selectedRegion.split(',')
			fig1 = output_plot_1(tumour_result, normal_result, tumour, TSA_list, outputPath, dpi)
			fig2 = output_plot_2(total_result, tumour, TSA_list, outputPath, dpi, withZero, selectedPeptide, selectedRegion)
			fig3 = output_plot_3(total_result, tumour, TSA_list, outputPath, dpi, withZero, selectedPeptide, selectedRegion)
			return [df1,df2,df3,fig1,fig2,fig3]
		else :
			fig1 = output_plot_1(tumour_result, normal_result, tumour, TSA_list, outputPath, dpi)
			fig2 = output_plot_2(total_result, tumour, TSA_list, outputPath, dpi, withZero, selectedPeptide, selectedRegion)
			fig3 = output_plot_3(total_result, tumour, TSA_list, outputPath, dpi, withZero, selectedPeptide, selectedRegion)
			return [df1,df2,df3,fig1,fig2,fig3]
	except ValueError as err:
		print(f"{err}", file = sys.stderr)
		return err 
	except Exception:
		print(f"An error occurred, please check your inputs", file = sys.stderr)
		return "An error occurred, please check your inputs"


def output_overall_csv(tumour_result, total_result, TSA_list, outputPath, withZero):
    data_rows = []  

    for peptide in total_result.keys():
        for region in total_result[peptide]:
            if region == 'total_count':
                continue
            for sample in total_result[peptide][region]['samples']:
                if int(withZero) == 0 and int(total_result[peptide][region]['samples'][sample]) == 0:
                    continue

                HERV_region = f"{total_result[peptide][region]['data']['chr']}:{total_result[peptide][region]['data']['start']}-{total_result[peptide][region]['data']['end']}"
                total_tumour_reads = tumour_result[peptide][region]['data']['total_count_region']
                total_reads = total_result[peptide][region]['data']['total_count_region']
                sample_total_reads = total_result[peptide][region]['samples'][sample]
                ORF = region
                TSA = TSA_list[int(peptide) - 1]
                strand = total_result[peptide][region]['data']['strand']
                
                # Get the padding sequence from the fasta file
                with open(f"{outputPath}/extracted_{peptide}_{region}_seq.fasta") as file:
                    for seq in SeqIO.parse(file, "fasta"):
                        padding_seq = str(seq.seq.upper())
                        if strand == '-':
                            padding_seq = str(seq.reverse_complement().seq.upper())

                row = {
                    "Sample": sample,
                    "HERV regions": HERV_region,
                    "Tumour reads": total_tumour_reads,
                    "Total reads": total_reads,
                    "Sample total reads": sample_total_reads,
                    "ORF": ORF,
                    "TSA": TSA,
                    "Strand": strand,
                    "Validation reading sequence": padding_seq
                }

                data_rows.append(row)

    # Create DataFrame from list of dictionaries
    df = pd.DataFrame(data_rows)

    df.to_csv(f"{outputPath}/Quantification_Summary_Across_All_Samples.csv", index=False)

    return df



def output_total_max_csv(tumour_result, total_result, TSA_list, outputPath):
	data_rows = []


	## get the region with maximum expression for each peptide
	max_region = {}

	for peptide in total_result:
		## i = the maximum number of count
		i = -1
	
		for region in total_result[peptide]:
			if region == "total_count":
				continue
			else:
				## if the looped region is larger than the saved i, update it
				if total_result[peptide][region]['data']['total_count_region'] > i:
					max_region[peptide] = \
						{
						"region_data" : total_result[peptide][region]["data"],
						"region_name" : region,
						}
					i = total_result[peptide][region]['data']['total_count_region']

		HERV_region = f"{max_region[peptide]['region_data']['chr']}:{max_region[peptide]['region_data']['start']}-{max_region[peptide]['region_data']['end']}"
		total_tumour_reads = tumour_result[peptide][max_region[peptide]["region_name"]]['data']['total_count_region']
		total_reads = total_result[peptide][max_region[peptide]["region_name"]]['data']['total_count_region']
		ORF = max_region[peptide]["region_name"]
		TSA = TSA_list[int(peptide) - 1]
		strand = max_region[peptide]["region_data"]['strand']
		## get the padding sequence from the fasta file
		with open(f"{outputPath}/extracted_{peptide}_{max_region[peptide]['region_name']}_seq.fasta") as file:
			for seq in SeqIO.parse(file, "fasta"):
				padding_seq = str(seq.seq.upper())
				if (max_region[peptide]['region_data']['strand'] == '-'):
					padding_seq = str(seq.reverse_complement().seq.upper())
		
		data_rows.append({
			"Peptide": peptide,
			"HERV regions": HERV_region,
			"Tumour reads": total_tumour_reads,
			"Total reads": total_reads,
			"ORF": ORF,
			"TSA": TSA,
			"Strand": strand,
			"Validation reading sequence": padding_seq
		})

	df = pd.DataFrame(data_rows)
	df.to_csv(f"{outputPath}/Maximal_HERV_Region_Counts_per_Query_Peptide_Across_All_Samples.csv", index=False)

	return df


def output_sample_max_csv(total_result, outputPath):
	
	max_region = {}
	sample_max = {}
	data_rows = []

	## get the sample ids from first region in first peptide
	total_result_iter = iter(total_result["1"])
	next(total_result_iter)
	first_region = next(total_result_iter)
	for sample in total_result["1"][first_region]['samples']:
		max_region[sample] = {}
		for peptide in total_result:
			max_region[sample][peptide] = \
				{
				'regions' : [],
				'counts' : []
				}

	## append all regions and corresponding count in two list for each sample
	for peptide in total_result:
		sample_max[peptide] = {}
		for region in total_result[peptide]:
			if (region == "total_count") :
				continue
			else :
				for sample in total_result[peptide][str(region)]['samples']:
					max_region[sample][str(peptide)]['regions'].append(region)
					max_region[sample][str(peptide)]['counts'].append(int(total_result[peptide][region]['samples'][sample]))

	## extract only the maximum
	for sample in max_region:
		for peptide in total_result:
			maximum = max(max_region[sample][peptide]['counts'])
			if (maximum == 0):
				sample_max[peptide][sample] = "no read"
			else:
				index = 0
				for count in max_region[sample][peptide]['counts']:
					if count == maximum:
						sample_max[peptide][sample] = max_region[sample][peptide]['regions'][index]
					index += 1


	for peptide in total_result:
		row = {
				"Peptide": peptide
			}
		for sample in sample_max[peptide]:
			row[sample] = sample_max[peptide][sample]
		data_rows.append(row)
	df = pd.DataFrame(data_rows)
	df.to_csv(f"{outputPath}/Maximal_HERV_Region_for_Each_Query_Peptide_by_Sample.csv", index=False)

	return df

def output_plot_1(tumour_result, normal_result, tumour, TSA_list, outputPath, dpi):
	
	i = 0

	## x_axis : tumour and normal
	grouping = ["Tumour", "Normal"]
	## y_axis : peptide
	peptides = []
	peptides_for_plot = []

	peptide_id = 1
	for peptide in TSA_list:
		if (str(peptide_id) in tumour_result):
			peptides.append(peptide_id)
			peptides_for_plot.append(TSA_list[int(peptide_id)-1])
		peptide_id+=1
	
	## get the total counts
	counts = []
	for peptide in peptides:
		for group in grouping:
			tumour_counts = tumour_result[str(peptide)]["total_count"]
			normal_counts = normal_result[str(peptide)]["total_count"]
		counts.append([tumour_counts, normal_counts])

	cell_height = 1.5

	fig, ax = plt.subplots(figsize=(2 * cell_height + 5, cell_height * len(peptides_for_plot) * 1.5 + 1))
	
	ax.set_title('Distribution of Total Counts in \nTumour and Normal Samples for Each Peptide', fontsize=20)
	ax.set_xlabel('Sample Classification', fontsize=20)
	ax.set_ylabel('Peptides', fontsize=20)
	
	im = ax.imshow(counts, cmap='Blues')
	divider = make_axes_locatable(ax)
	cax = divider.append_axes('right', size='5%', pad=0.1)
	cbar = ax.figure.colorbar(im, cax=cax, cmap='Blues', shrink = 0.5)
	cbar.set_label('counts', rotation = -90,labelpad=10 , fontsize=20)
	
	ax.set_yticks(np.arange(len(peptides)) , labels=peptides_for_plot, minor = False)
	ax.set_xticks(np.arange(len(grouping)), labels=grouping)

	plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
	plt.setp(ax.get_yticklabels(), rotation=45, va='center', ha="right", rotation_mode="anchor")

	plt.savefig(f"{outputPath}/Distribution_of_Total_Counts_in_Tumour_and_Normal_Samples_for_Each_Peptide.png", dpi=dpi)

	return fig



def output_plot_2(total_result, tumour, TSA_list, outputPath, dpi, withZero, selectedPeptide, selectedRegion):
	figures = []
	i = 0
	peptide_len = len(total_result.keys())

	## get the list of samples and regions
	## regions save information, regions_for_plot save simplified name
	regions = {}
	regions_for_plot = {}

	## x axis : tumour or normal sample
	grouping = ["Tumour", "Normal"]

	for peptide in total_result:
		regions_for_plot[peptide] = []
		regions[peptide] = {}
		for region in total_result[peptide]:
			if (region == 'total_count'):
				continue
			else:
				regions[peptide]["_".join(region.split("_")[0:3])] = region
				regions_for_plot[peptide].append("_".join(region.split("_")[0:3]))
	

	## if user only need specific peptide and region
	if (selectedPeptide != None):
		regions_for_plot = {
			selectedPeptide : regions_for_plot[selectedPeptide]
		}
		regions = {
			selectedPeptide : regions[selectedPeptide]
		}
		if (selectedRegion != None):
			regions_for_plot = {
				selectedPeptide : selectedRegion
			}
			r_list = list(regions[selectedPeptide].keys())
			for r in r_list:
				if r not in selectedRegion:
					regions[selectedPeptide].pop(r)

	
	## counts : peptide (list) => each region (lists) => group's total count (integer)
	counts = {}
	for peptide in regions:
		counts[peptide] = []
		i = 0
		for region in regions[peptide]:
			tumour_count = 0
			normal_count = 0
			
			## get the corresponding count and save in 2d array	
			for sample in tumour:
				tumour_count += int(total_result[peptide][regions[peptide][region]]['samples'][f"{sample}T"])
				normal_count += int(total_result[peptide][regions[peptide][region]]['samples'][f"{sample}N"])
			
			## filter zero if withZero == 0
			if (withZero == 0):
				if (tumour_count == 0 and normal_count == 0):
					if (region in regions_for_plot[peptide]):
						regions_for_plot[peptide].remove(region)
					continue

			counts[peptide].append([])
			counts[peptide][i].append(tumour_count)
			counts[peptide][i].append(normal_count)
			i += 1

	

	for peptide in regions_for_plot:
		cell_height = 1.5	
		fig, ax = plt.subplots(figsize=(10 ,len(regions_for_plot[peptide])*cell_height+3))
		
		ax.set_title(f"Distribution of Total Counts in Tumour and Normal\nSamples Across HERV Regions of Peptide {peptide}:{TSA_list[int(peptide)-1]}", fontsize=20)
		ax.set_xlabel('Sample Classification', fontsize=20)
		ax.set_ylabel('Regions', fontsize=20)

		im = ax.imshow(counts[peptide], cmap='Blues')

		divider = make_axes_locatable(ax)
		cax = divider.append_axes('right', size='5%', pad=0.1)
		cbar = ax.figure.colorbar(im, cax=cax, cmap='Blues', orientation = 'vertical')
		cbar.set_label('counts', labelpad=10 , fontsize=20)
		
		ax.set_yticks(np.arange(len(regions_for_plot[peptide])), labels=regions_for_plot[peptide])
		ax.set_xticks(np.arange(len(grouping)), labels=grouping)

		plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
		plt.setp(ax.get_yticklabels(), rotation=45, ha="right", rotation_mode="anchor")

		plt.savefig(f"{outputPath}/Distribution_of_Total_Counts_in_Tumour_and_Normal_Samples_Across_HERV_Regions_of_Peptide_{peptide}.png", dpi=dpi)

		figures.append(fig)

	return figures



def output_plot_3(total_result, tumour, TSA_list, outputPath, dpi, withZero, selectedPeptide, selectedRegion):
	figures = []
	i = 0
	peptide_len = len(total_result.keys())

	## get the list of samples and regions
	## regions save information, regions_for_plot save simplified name
	regions = {}
	regions_for_plot = {}
	
	## samples save sample as a pair, samples_for_plot save seperately
	samples = []
	samples_for_plot = []
	
	for peptide in total_result:
		regions_for_plot[peptide] = []
		regions[peptide] = {}
		for region in total_result[peptide]:
			if (region == 'total_count'):
				continue
			else:
				regions[peptide]["_".join(region.split("_")[0:3])] = region
				regions_for_plot[peptide].append("_".join(region.split("_")[0:3]))
	
	for sample in tumour:
		samples.append([f"{sample}T", f"{sample}N"])
		samples_for_plot.append(f"{sample}T")
		samples_for_plot.append(f"{sample}N")

	## if user only need specific region and peptide
	if (selectedPeptide != None):
		regions_for_plot = {
			selectedPeptide : regions_for_plot[selectedPeptide]
		}
		regions = {
			selectedPeptide : regions[selectedPeptide]
		}
		if (selectedRegion != None):
			regions_for_plot = {
				selectedPeptide : selectedRegion
			}
			r_list = list(regions[selectedPeptide].keys())
			for r in r_list:
				if r not in selectedRegion:
					regions[selectedPeptide].pop(r)

	## counts : peptide (list) => each region (lists) => sample's count (integer)
	counts = {}
	for peptide in regions:
		counts[peptide] = []
		i = 0
		for region in regions[peptide]:
			counts[peptide].append([])

			## get the corresponding count and save in 2d array	
			for sample in samples_for_plot:
				sample_count = total_result[peptide][regions[peptide][region]]['samples'][sample]
				counts[peptide][i].append(int(sample_count))

			## remove zero if withZero == 0
			if (withZero == 0):
				if (all(item == 0 for item in counts[peptide][i])):
					if (region in regions_for_plot[peptide]):
						regions_for_plot[peptide].remove(region)
						counts[peptide].remove(counts[peptide][i])
					continue
			i += 1


	for peptide in regions_for_plot:

		fig, ax = plt.subplots(figsize=(len(samples_for_plot) * 0.75 + 2.5, len(regions_for_plot[peptide]) * 0.75+3))
		ax.set_title(f"Distribution of Total Counts Across HERV Regions \nin Each Sample for Peptide {peptide}:{TSA_list[int(peptide)-1]}", fontsize=20)
		ax.set_xlabel('Samples', fontsize=20)
		ax.set_ylabel('Regions', fontsize=20)

		im = ax.imshow(counts[peptide], cmap='Blues')

		divider = make_axes_locatable(ax)
		cax = divider.append_axes('right', size='5%', pad=0.1)
		cbar = ax.figure.colorbar(im, cax=cax, cmap='Blues', orientation = 'vertical')
		cbar.set_label('counts', labelpad=10, fontsize=20)
		
		ax.set_yticks(np.arange(len(regions_for_plot[peptide])), labels=regions_for_plot[peptide])
		ax.set_xticks(np.arange(len(samples_for_plot)), labels=samples_for_plot)

		plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
		plt.setp(ax.get_yticklabels(), rotation=45, ha="right", rotation_mode="anchor")

		plt.savefig(f"{outputPath}/Distribution_of_Total_Counts_Across_HERV_Regions_in_Each_Sample_for_Peptide_{peptide}.png", dpi=dpi)
		figures.append(fig)

	return figures


def getTSAlist(inputPeptide):
	## get the input peptide sequence
	with open(inputPeptide, 'r') as file:
		TSA_list = []
		for line in file:
			if ">" in line:
				continue
			line = line.replace("\n","")
			TSA_list.append(line)
	return TSA_list


def extractResult(sampleResultDir, TSA_list):

	with open(sampleResultDir, 'r') as file:
		result_list = file.readlines()

	result_list = [line.strip().split(',') for line in result_list]
	sample_result_dir = []

	for result in result_list:
		if (not "summary" in result[2]):
			sample_result_dir.append(result[2])

	sample_result_dir.pop(0)

	result_summary = {}
	HERV_peptides = []
	i = 1
	
	for peptide in TSA_list:
		result_summary[str(i)] = {"total_count": 0}
		i += 1
	
	with open(sample_result_dir[0], 'r') as file:
		## skip the first two lines
		next(file)
		next(file)

		for line in file:
			fields = line.strip().split('\t')
			peptide_seq = fields[0].split('_')[1]
			region = fields[0]
			sample = sample_result_dir[0].split("/")[-1].split("_")[0]
			result_summary[peptide_seq][region] = \
				{
				"data" : 
					{
					"chr" : fields[1],
					"start" : fields[2],
					"end" : fields[3],
					"strand" : fields[4],
					"length" : fields[5],
					"total_count_region" : int(fields[6])
					},
				"samples" :
					{
					sample : fields[6]
					}
				}
			result_summary[peptide_seq]["total_count"] += int(fields[6])

			if (peptide_seq not in HERV_peptides):
				HERV_peptides.append(peptide_seq)
			
	sample_result_dir.pop(0)

	for result_dir in sample_result_dir:

		with open(result_dir, 'r') as file:
			next(file)
			next(file)

			for line in file:
				fields = line.strip().split('\t')
				peptide_seq = fields[0].split('_')[1]
				region = fields[0]
				sample = result_dir.split("/")[-1].split("_")[0]
				result_summary[peptide_seq][region]["data"]["total_count_region"] += int(fields[6])
				result_summary[peptide_seq][region]["samples"][sample] = fields[6]
				result_summary[peptide_seq]["total_count"] += int(fields[6])

	for peptide in list(result_summary):
		if peptide not in HERV_peptides:
			result_summary.pop(peptide)

	return result_summary


def combine_tumour_normal_results(tumour_result, normal_result):

	total_result = {}
	for peptide in tumour_result:
		total_result[peptide] = {}
		for region in tumour_result[peptide]:
			if region == "total_count":
				total_result[peptide]["total_count"] = tumour_result[peptide]["total_count"] + normal_result[peptide]["total_count"]
			else:
				total_result[peptide][region] = copy.deepcopy(tumour_result[peptide][region])
				total_result[peptide][region]["data"]["total_count_region"] = tumour_result[peptide][region]["data"]["total_count_region"] + normal_result[peptide][region]["data"]["total_count_region"]
				for sample_name in tumour_result[peptide][region]["samples"]:
					total_result[peptide][region]["samples"] = {**normal_result[peptide][region]["samples"], **tumour_result[peptide][region]["samples"]}
	
	return total_result

def get_padding(total_result, outputPath):

	## produce the bed file
	for peptide in total_result.keys():
		for region in total_result[peptide]:

			if (region == 'total_count'):
				continue
			
			bed_start = int(total_result[peptide][region]['data']['start']) - 60
			bed_end = int(total_result[peptide][region]['data']['end']) + 60

			with open(f"{outputPath}/{peptide}_{region}.bed", 'w') as bed_file:
				bed_line = (
					f"{total_result[peptide][region]['data']['chr']}\t"
					f"{bed_start}\t"
					f"{bed_end}\t"
					f"{peptide}_regions.fasta"
					)
				bed_file.write(bed_line + "\n")

			bedCommand = f"bedtools getfasta -name -fi ./data/hg19/{total_result[peptide][region]['data']['chr']}.fa -bed {outputPath}/{peptide}_{region}.bed -fo {outputPath}/extracted_{peptide}_{region}_seq.fasta"
			try:
				subprocess.run(bedCommand, shell = True, check = True)
			except subprocess.CalledProcessError as err:
				return err

	return None










  
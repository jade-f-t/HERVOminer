import subprocess
import os 
import sys
import json 
import pysam
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

def func1_blastORF(inputPeptide, outputPath, inputVCF):
	"""
	Find regions of similarity between query peptides and the open reading frames of HERV regions.
	Inputs :
		inputPeptide : path to the fasta file containing the peptide sequences
		outputPath : path of the output file
	Output : blastp_output.txt
	"""
	if inputVCF is None:
		blastDBdir = "./data/HervOrfBlastpDB/HervOrfBlastpDB"
	else:
		[blastDBdir, makeDBerr] = individualVCF(inputVCF, outputPath)
	
	blastError = blastp(inputPeptide, outputPath, blastDBdir)
	
	if (blastError is not None):
		print(f"Error: please check your inputs: {blastError}", file = sys.stderr)
		raise ValueError(f"Error: please check your inputs: {blastError}")
		sys.exit(1)
	
	return None

def blastp(inputPeptide, outputPath, blastDBdir):
	blastpCommand = [
		'blastp', 
		'-query', inputPeptide, '-db', blastDBdir,
		'-word_size', '3', '-gapopen', '9', '-gapextend', '1', '-matrix', 'PAM30', '-threshold', '16',
		'-comp_based_stats', '0', '-window_size', '15', '-num_threads', '10', '-evalue', '0.05',
		'-outfmt', "'6 qacc sacc pident qseq length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore'",
		'-out', os.path.join(outputPath, 'blastp_output.txt')
	]

	try:
		subprocess.run(' '.join(blastpCommand),shell = True, check = True)
		return None
	except subprocess.CalledProcessError as err:
		return err


def individualVCF(inputVCF, outputPath):
	
	HERV_ORF_dict_dir = "/home/ftw0710/HERVOminer_server/data/HERV_ORF_dict"
	HERV_ORF_bed_dir = "/home/ftw0710/HERVOminer_server/data/ORF_region.bed"
	HERV_ref_dir = "/home/ftw0710/HERVOminer_server/data/1_separatedSeq"
	mergedORF_dir = "/home/ftw0710/HERVOminer_server/data/3_mergeORF.fasta"
	vcf_dir = inputVCF
	vcf_bed_dir = f"{outputPath}/vcf_region.bed"
	output_dir = outputPath

	vcf_info = {}
	orf_id_list = []
	orf_id_vcf_id = {}

	try:
		# S1 : read ORF database
		with open(HERV_ORF_dict_dir, "r") as file:
			HERV_ORF_dict = json.load(file)

		# S2 : get vcf information, produce bed file for vcf 
		with open(vcf_bed_dir, "w") as file:
			with pysam.VariantFile(vcf_dir) as vcf_in:
				id = 0
				for record in vcf_in:
					chrom = record.chrom
					pos = record.pos
					ref = record.ref
					alts = record.alts[0]

					if len(ref) == len(alts):
						type = "snv"
						end = pos+1
					else:
						continue
					
					vcf_info[id] = {
						"start_position":pos,
						"end_position":end,
						"type":type,
						"ref":ref,
						"alt":alts
					}
					
					line = f"{chrom}\t{pos - 1}\t{end - 1}\t{id}\n"
					file.write(line)
					id+=1
		# S3 : check if variant intersect orf database
		bedCommand = f"bedtools intersect -wa -wb -a {HERV_ORF_bed_dir} -b {vcf_bed_dir} > {output_dir}/ORF_vcf_overlap.tsv"

		subprocess.run(bedCommand, shell = True, check = True)

		# S4 : change the ORF sequence in mergedORF (database)
		ORF_vcf_headers = ["chrom", "orf_start", "orf_end", "orf_id", "blank", "strand", "chromo_2", "vcf_start", "vcf_end", "vcf_id"]
		ORF_vcf_overlap_df = pd.read_csv(f"{output_dir}/ORF_vcf_overlap.tsv", sep="\t", header=None, names=ORF_vcf_headers)

		## S4.1 : get the rna sequence of the orf related
		for index, row in ORF_vcf_overlap_df.iterrows():

			# get the HERV sequence
			file_id = "_".join(row['orf_id'].split("_")[:3])
			HERV_seq_file = f"{HERV_ref_dir}/{file_id}.fasta"
			record = next(SeqIO.parse(HERV_seq_file, "fasta"))
			description = record.description.split('\t')
			strand = description[1]

			sequence = str(record.seq)
			if strand == '-':
				sequence = str(record.seq.reverse_complement())

			## S4.2 : change the HERV sequence 
			HERV_region = row['orf_id'].split("_")[2]
			HERV_start = int(HERV_region.split("-")[0])
			HERV_end = int(HERV_region.split("-")[1])

			vcf_start = row['vcf_start'] + 1
			vcf_id = row['vcf_id']
			chrom = row['chrom']

			position = int(vcf_start) - int(HERV_start)
			new_sequence_list = list(sequence)
			
			new_sequence_list[position] = vcf_info[vcf_id]['alt']

			new_sequence = str("".join(new_sequence_list)).upper()
			original_sequence = sequence
			if strand == "-":
				new_sequence = str(Seq("".join(new_sequence_list)).reverse_complement()).upper()
				original_sequence = str(Seq(sequence).reverse_complement()).upper()

			## S4.3 : find which part of the sequence responsible for the ORF
			orf_start = int(row['orf_start']) + 1
			orf_end = int(row['orf_end'])

			if strand == "-":
				start = HERV_end - orf_end # should be +1 but consider start from 0
				end = orf_end - orf_start + start # position of orf end and orf start is according to the forward strand so end > start
			else:
				start = orf_start - HERV_start
				end = orf_end - orf_start + start 
			# record the orf (new) sequence in nucleotide and the orf id
			if 'orf' not in vcf_info[vcf_id].keys():
				vcf_info[vcf_id]['orf'] = {
					row['orf_id']:new_sequence[start:end+1]
				}
				vcf_info[vcf_id]['orf_region'] = {
					row['orf_id']:[orf_start, orf_end]
				}
				vcf_info[vcf_id]['orf_strand'] = {
					row['orf_id']:strand
				}
				vcf_info[vcf_id]['orf_original_seq'] = {
					row['orf_id']:original_sequence[start:end+1]
				}
			else:
				vcf_info[vcf_id]['orf'][row['orf_id']] = new_sequence[start:end+1]
				vcf_info[vcf_id]['orf_region'][row['orf_id']] = [orf_start, orf_end]
				vcf_info[vcf_id]['orf_strand'][row['orf_id']] = strand
				vcf_info[vcf_id]['orf_original_seq'][row['orf_id']] = original_sequence[start:end+1]
			orf_id_list.append(row['orf_id'])
			if row['orf_id'] not in orf_id_vcf_id.keys():
				orf_id_vcf_id[row['orf_id']] = [vcf_id]
			else:
				orf_id_vcf_id[row['orf_id']].append(vcf_id)
		## S4.4 : create a updated merge ORF fasta according to the variants (indel not yet) 
		with open(f"{output_dir}/merged_orf_variant.fasta", 'w') as outfile:
			for record in SeqIO.parse(mergedORF_dir, "fasta"):
				description = record.description
				sequence_id = str(record.id)
				orf_sequence = str(record.seq)
				

				outfile.write(f">{description}\n")
				
				if sequence_id not in orf_id_list:
					aa_sequence = str(orf_sequence)
				else:
					vcf_ids = orf_id_vcf_id[sequence_id]
					for vcf_id in vcf_ids:
						update_sequence = list(vcf_info[vcf_id]['orf_original_seq'][sequence_id])
						original_sequence = list(vcf_info[vcf_id]['orf_original_seq'][sequence_id])
						for vcf_id in orf_id_vcf_id[sequence_id]:
							sequence_variant = list(vcf_info[vcf_id]['orf'][sequence_id])
						
							# if vcf_info[vcf_id]['type'] == 'snv':
							if vcf_info[vcf_id]['orf_strand'][sequence_id] == '-':
								pos = vcf_info[vcf_id]['orf_region'][sequence_id][1] - vcf_info[vcf_id]['start_position']
								update_sequence[pos] = str(Seq(sequence_variant[pos]).reverse_complement()).upper()

							else:
								pos = vcf_info[vcf_id]['start_position'] - vcf_info[vcf_id]['orf_region'][sequence_id][0]
								update_sequence[pos] = sequence_variant[pos]
							
						final_sequence = "".join(update_sequence)
						ORF_frame = int(sequence_id.split('_')[-2])

						if ORF_frame > 3:
							final_sequence = str(Seq(final_sequence).reverse_complement())

							aa_sequence = str(Seq(final_sequence).translate())
						if ORF_frame <= 3:
							aa_sequence = str(Seq(final_sequence).translate())
						
						original = HERV_ORF_dict[sequence_id]["ORF"]
						
						# S5 : adjust the ORF sequence in HERV_ORF_dict
						if original != aa_sequence:
							HERV_ORF_dict[sequence_id]["ORF"] = aa_sequence
				outfile.write(aa_sequence + "\n")

		with open(f"{output_dir}/HERV_ORF_dict_variant", 'w') as outfile:
				json.dump(
					HERV_ORF_dict,  
					outfile,        
					indent=4         
				)

		# S6 : build blastp database by the command: makeblastdb -in 3_mergeORF.fasta -input_type fasta -dbtype prot -parse_seqids -title "HERV ORF blastp Database" -out HervOrfBlastpDB
		makeblastdb_command = f"makeblastdb -in {output_dir}/merged_orf_variant.fasta -input_type fasta -dbtype prot -parse_seqids -title 'HERV ORF variant blastp Database' -out {output_dir}/HervOrfBlastpDB_variant"
		subprocess.run(makeblastdb_command, shell = True, check = True)

		makeDBerr = None
	except subprocess.CalledProcessError as err:
		makeDBerr = error
	return [f"{output_dir}/HervOrfBlastpDB_variant",makeDBerr]
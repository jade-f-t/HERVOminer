import subprocess
import os 
import sys

def func1_blastORF(inputPeptide, outputPath):
			
	blastError = blastp(inputPeptide, outputPath)
	if (blastError is not None):
		print('Error: please check your inputs', file = sys.stderr)
		return blastError
		sys.exit(1)

	return None

def blastp(inputPeptide, outputPath):
	blastpCommand = \
		f"blastp -query {inputPeptide} \
		-db ./data/HervOrfBlastpDB/HervOrfBlastpDB \
		-word_size 3 -gapopen 9 -gapextend 1 -matrix PAM30 -threshold 16 \
		-comp_based_stats 0 -window_size 15 -num_threads 48 -evalue 0.05 \
		-outfmt '6 qacc sacc pident qseq length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore'\
		-out {outputPath}/blastp_output.txt"

	try:
		subprocess.run(blastpCommand, shell = True, check = True)
		return None
	except subprocess.CalledProcessError as err:
		return err




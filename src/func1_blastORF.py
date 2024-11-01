import subprocess
import os 
import sys

def func1_blastORF(inputPeptide, outputPath):
	"""
	Find regions of similarity between query peptides and the open reading frames of HERV regions.
	Inputs :
		inputPeptide : path to the fasta file containing the peptide sequences
		outputPath : path of the output file
	Output : blastp_output.txt
	"""
	blastError = blastp(inputPeptide, outputPath)
	
	if (blastError is not None):
		print(f"Error: please check your inputs: {blastError}", file = sys.stderr)
		raise ValueError(f"Error: please check your inputs: {blastError}")
		sys.exit(1)
	
	return None

def blastp(inputPeptide, outputPath):
	blastpCommand = [
		'blastp', 
		'-query', inputPeptide, '-db', './data/HervOrfBlastpDB/HervOrfBlastpDB',
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




import argparse
import os 

def build_parser():
	parser = argparse.ArgumentParser(prog = "HERVtool", usage = "... to be added", description = "... to be added")
	subcmd = parser.add_subparsers(dest='subcmd', help='subcommands', metavar='SUBCOMMAND')
	subcmd.required = True

	# function 1 blastORF (func1_blastORF.py) : 
	blastORF_parser = subcmd.add_parser('blastORF', help = 'blastp the input file to the orf database')
	blastORF_parser.add_argument("-p", "--peptide", required=True, help = "path to the input peptide file", dest = "inputPeptide")
	blastORF_parser.add_argument("-o", "--output", default = os.getcwd(), help = "path to the output file", dest = "outputPath")

	# function 2 summarize and annotation (func2_summarize_annotation.py) :
	summarize_annotation_parser = subcmd.add_parser('summarizeAnnotate', help = 'get the corresponding HERV region compare to the blastp output result get from function 1 and create annotation file in gtf format')
	summarize_annotation_parser.add_argument("-i", "--input_blastpOutput", required=True, help = "path to the input blastp output file get from function 1", dest = "blastpOutput")
	summarize_annotation_parser.add_argument("-o", "--output", default = os.getcwd(), help = "path to the output file", dest = "outputPath")

	# # function 3 createGTF (func3_createGTF.py) :
	# createGTF_parser = subcmd.add_parser('createGTF', help = 'create annotation file in gtf format')
	# createGTF_parser.add_argument("-i", "--input_HERVregion", required = True, help = "path to the HERV region information json file get from function 2", dest = "HERVregion")
	# createGTF_parser.add_argument("-o", "--output", default = os.getcwd(), help = "path to the output file", dest = "outputPath")

	# function 3 quantification (func3_quantification.py) :
	quantification_parser = subcmd.add_parser('quantification', help = 'use featureCount count to quantify the HERV region')
	quantification_parser.add_argument("-i", "--input_csv_file", required=True, help = "path to the input csv file with all of the bam file path", dest = "inputCsvFile")
	quantification_parser.add_argument("-a", "--annotation", required = True, help = "path to the annotation file created by function 2", dest = "annotationFile")
	quantification_parser.add_argument("-o", "--output", default=os.getcwd(), help = "path to the output file", dest = "outputPath")
	quantification_parser.add_argument("-t", "--threads", required = True, help = "no of threads to use for one sample", dest = "threadNo" )
	quantification_parser.add_argument("-n", "--parallel_task", required = True, help = "number of sample to be analyzed in parallel computing", dest = "parallelTask")

	# function 4 outputResult (func4_outputResult.py) :
	outputResult_parser = subcmd.add_parser('outputResult', help = 'output the analyzed result in csv format')
	outputResult_parser.add_argument("-p", "--peptide", required=True, help = "path to the input peptide file", dest = "inputPeptide")
	outputResult_parser.add_argument("-i", "--input_csv_file", required=True, help = "path to the input csv file with all of the bam file path", dest = "inputCsvFile")
	outputResult_parser.add_argument("-T", "--tumour_quantification_result", required=True, help = "path to the resultDirectories_T.csv file", dest = "tumourResultDirectories")
	outputResult_parser.add_argument("-N", "--normal_quantification_result", required=True, help = "path to the resultDirectories_N.csv file", dest = "normalResultDirectories")
	outputResult_parser.add_argument("-z", "--with_zero", default = 0, help = "include region with zero count or not, 1 : with zero, 0: without zero", dest = "withZero")
	outputResult_parser.add_argument("-o", "--output", default = os.getcwd(), help = "path to the output file", dest = "outputPath")
	outputResult_parser.add_argument("-d", "--dpi", default = 100, help = "set the dpi of the figures", dest = "dpi")


	# HERVtool(args.inputPeptide, args.inputCsvFile, args.outputPath, args.thread, args.group)
	HERVOminer_parser = subcmd.add_parser('HERVOminer', help = 'finish the whole process')
	HERVOminer_parser.add_argument("-p", "--peptide", required=True, help = "path to the input peptide file", dest = "inputPeptide")
	HERVOminer_parser.add_argument("-i", "--input_csv_file", required=True, help = "path to the input csv file with all of the bam file path", dest = "inputCsvFile")
	HERVOminer_parser.add_argument("-o", "--output", default = os.getcwd(), help = "path to the output file", dest = "outputPath")
	HERVOminer_parser.add_argument("-t", "--threads", required = True, help = "no of threads to use for one sample", dest = "threadNo" )
	HERVOminer_parser.add_argument("-n", "--parallel_task", required = True, help = "number of sample to be analyzed in parallel computing", dest = "parallelTask")
	HERVOminer_parser.add_argument("-z", "--with_zero", default = 0, help = "include region with zero count or not, 1 : with zero, 0: without zero", dest = "withZero")
	HERVOminer_parser.add_argument("-d", "--dpi", default = 100, help = "set the dpi of the figures", dest = "dpi")


	return parser 




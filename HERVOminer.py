#!/usr/bin/env python

import argparse
import sys

from build_parser import build_parser

from ./src/func1_blastORF import func1_blastORF
from ./src/func2_summarize_annotation import func2_summarize_annotation
from ./src/func3_quantification import func3_quantification
from ./src/func4_outputResult import func4_outputResult

if __name__ == "__main__":
	parser = build_parser()
	args = parser.parse_args()

	if args.subcmd == 'blastORF':
		returnValue = func1_blastORF(args.inputPeptide, args.outputPath)

	elif args.subcmd == 'summarizeAnnotate':
		returnValue = func2_summarize_annotation(args.blastpOutput, args.outputPath)

	# elif args.subcmd == 'createGTF':
	# 	returnValue = func3_createGTF(args.HERVregion, args.outputPath)

	elif args.subcmd == 'quantification':
		returnValue = func3_quantification(args.inputCsvFile, args.annotationFile, args.outputPath,args.threadNo, args.parallelTask)

	elif args.subcmd == 'outputResult':
		returnValue = func4_outputResult(args.inputPeptide, args.inputCsvFile, args.tumourResultDirectories, args.normalResultDirectories, args.outputPath, args.withZero, args.dpi)

	elif args.subcmd == 'HERVOminer':
		blastpOutput = f"{args.outputPath}/blastp_output.txt"
		annotationFile = f"{args.outputPath}/quantification.gtf"
		tumourResultDirectories = f"{args.outputPath}/resultDirectories_T.csv"
		normalResultDirectories = f"{args.outputPath}/resultDirectories_N.csv"
		func1_blastORF(args.inputPeptide, args.outputPath)
		func2_summarize_annotation(blastpOutput, args.outputPath)
		func3_quantification(args.inputCsvFile, annotationFile, args.outputPath,args.threadNo, args.parallelTask)
		func4_outputResult(args.inputPeptide, args.inputCsvFile, tumourResultDirectories, normalResultDirectories, args.outputPath, args.withZero, args.dpi)

		returnValue = None

	sys.exit(returnValue)



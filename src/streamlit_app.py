import streamlit as st
import traceback
from func1_blastORF import func1_blastORF
from func2_summarize_annotation import func2_summarize_annotation
from func3_quantification import func3_quantification
from func4_outputResult import func4_outputResult

## test data
p = '/work1783/ftw0710/0117test/HERV_fasta.fa'
i = '/work1786/user_home_2/ftw0710/HERVpipeline/commandLine/input.csv'
o =  '/work1783/ftw0710/package_test'

def dataChart(s4results):
	df1, df2, df3, fig1, fig2, fig3 = st.tabs(["Table 1", "Table 2", "Table 3", "Plot 1", "Plot 2", "Plot 3"])
	with df1:
		st.header("Quantification Summary Across All Samples")
		st.dataframe(s4results[0])
	with df2:
		st.header("Maximal HERV Region Counts per Query Peptide Across All Samples")
		st.dataframe(s4results[1])
	with df3:
		st.header("Maximal HERV Region for Each Query Peptide by Sample")
		st.dataframe(s4results[2])
	with fig1:
		st.pyplot(s4results[3])
	with fig2:
		for i in range(len(s4results[4])):
			st.pyplot(s4results[4][i])
	with fig3:
		for i in range(len(s4results[5])):
			st.pyplot(s4results[5][i])

st.title("HERVOminer")

StepByStep, Integrated = st.tabs(["Step-by-Step Process", "Integrated Process"]) 

with StepByStep:
	with st.form("blastORF"):
		
		st.write("step 1 : Protein BLAST Similarity Analysis (BLAST)")
		inputPeptide = st.text_input('path to the input peptide file', p)
		p = inputPeptide
		outputPath = st.text_input('output path', o)
		o = outputPath
		
		submitted = st.form_submit_button("Submit")
		if (submitted):
			try:
				with st.spinner('Wait for it...'):
					err = func1_blastORF(inputPeptide, outputPath)
					if (err != None):
						raise ValueError(err)
				st.success('Step 1 completed!')
			except ValueError as err:
				st.error(f"An error occurred, please check your inputs:\n{err}")
			except Exception:
				err = traceback.format_exc()
				st.error(f"An error occurred, please check your inputs:\n{err}")

	st.divider()

	with st.form("summarizeAnnotate"):
		
		st.write("step 2 : BLAST results summarisation, Annotation file generation")
		blastpOutput = st.text_input("path to the input blastp output file get from Step 1", f"{o}/blastp_output.txt")
		outputPath = st.text_input('output path', o)
		o = outputPath
		
		submitted = st.form_submit_button("Submit")
		if (submitted):
			try:
				with st.spinner('Wait for it...'):
					func2_summarize_annotation(blastpOutput, outputPath)
				st.success('Step 2 completed!')
			except Exception:
				err = traceback.format_exc()
				st.error(f"An error occurred, please check your inputs:\n{err}")

	st.divider()

	with st.form("quantification"):
		st.write("step 3 : Quantification (featureCounts)")
		inputCsvFile = st.text_input("path to the input csv file with all of the bam file path", i)
		i = inputCsvFile
		annotationFile = st.text_input("path to the annotation file created by Step 2", f"{o}/quantification.gtf")
		outputPath = st.text_input('output path', o)
		o = outputPath
		threadNo = st.number_input('no of threads to use for one sample', min_value = 1, value = 1)
		parallelTask = st.number_input("number of sample to be analyzed in parallel computing", min_value = 1, value = 1)
		
		submitted = st.form_submit_button("Submit")
		if (submitted):
			try:
				with st.spinner('Wait for it...'):
					func3_quantification(inputCsvFile, annotationFile, outputPath, threadNo, parallelTask)
				st.success('Step 3 completed!')
			except Exception:
				err = traceback.format_exc()
				st.error(f"An error occurred, please check your inputs:\n{err}")

	st.divider()

	with st.form("outputResult"):
		st.write("step 4 : Data summarisation and visualization")
		inputPeptide = st.text_input('path to the input peptide file', p)
		p = inputPeptide
		inputCsvFile = st.text_input("path to the input csv file with all of the bam file path", i)
		i = inputCsvFile
		tumourResultDirectories = st.text_input("path to the resultDirectories_T.csv file", f"{o}/resultDirectories_T.csv")
		normalResultDirectories = st.text_input("path to the resultDirectories_N.csv file", f"{o}/resultDirectories_N.csv")
		outputPath = st.text_input('output path', o)
		o = outputPath
		zeroCountException = st.toggle('Include region with zero count')
		if (zeroCountException):
			withZero = 1
		else:
			withZero = 0
		dpi = st.number_input('dpi of the figures', min_value = 0, value = 100)
		selectedPeptide = st.text_input("selected peptide to generate respective plots, please input the id of the peptide (optional)", None)
		selectedRegion = st.text_input("selected regions to appear in the plots, format : <region_id>,<region_id>,... eg. 2_1_3166,10_1_3514,15_1_2630 (optional)", None)

		submitted = st.form_submit_button("Submit")
		if (submitted):
			try:
				with st.spinner('Wait for it...'):
					s4results = func4_outputResult(inputPeptide, inputCsvFile, tumourResultDirectories, normalResultDirectories, outputPath, withZero, dpi, selectedPeptide, selectedRegion)
					if (type(s4results) != list):
						raise ValueError(s4results)
				st.success('Step 4 completed!')
				dataChart(s4results)
			except ValueError as err:
				st.error(f"An error occurred, please check your inputs:\n{err}")
			except Exception:
				err = traceback.format_exc()
				st.error(f"An error occurred, please check your inputs:\n{err}")

with Integrated:
	with st.form("Integrated Process"):
		st.header("Integrated Process")
		inputPeptide = st.text_input('path to the input peptide file', p)
		p = inputPeptide
		inputCsvFile = st.text_input("path to the input csv file with all of the bam file path", i)
		i = inputCsvFile
		outputPath = st.text_input('output path', o)
		o = outputPath
		threadNo = st.number_input('no of threads to use for one sample', min_value = 1, value = 1)
		parallelTask = st.number_input("number of sample to be analyzed in parallel computing", min_value = 1, value = 1)
		zeroCountException = st.toggle('Include region with zero count')
		if (zeroCountException):
			withZero = 1
		else:
			withZero = 0
		dpi = st.number_input('dpi of the figures', min_value = 0, value = 100)
		
		submitted = st.form_submit_button("Submit")
		if (submitted):
			try:
				with st.status("HERVOminer started...", expanded = True) as status:
					err = func1_blastORF(inputPeptide, outputPath)
					if (err != None):
						raise ValueError(err)
					st.write('Step 1 : Protein BLAST Similarity Analysis completed!')
					func2_summarize_annotation(blastpOutput, outputPath)
					st.write('Step 2 : BLAST results summarisation, Annotation file generation completed!')
					func3_quantification(inputCsvFile, annotationFile, outputPath, threadNo, parallelTask)
					st.write('Step 3 : Quantification completed')
					s4results = func4_outputResult(inputPeptide, inputCsvFile, tumourResultDirectories, normalResultDirectories, outputPath, withZero, dpi, selectedPeptide, selectedRegion)
					if (type(s4results) != list):
						raise ValueError(s4results)
					st.write('Step 4 : Data summarisation and visualization completed!')
					status.update(label="HERVOminer completed!", state="complete", expanded=False)
					dataChart(s4results)
			except ValueError as err:
				st.error(f"An error occurred, please check your inputs:\n{err}")
			except Exception:
				err = traceback.format_exc()
				st.error(f"An error occurred, please check your inputs:\n{err}")





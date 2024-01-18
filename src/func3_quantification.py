import subprocess
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

def func3_quantification(inputCsvFile, annotationFile, outputPath, threadNo, parallelTask):
	
	## get bam directories from csv input
	inputCsvResult, inputCsvError = inputFileHandle(inputCsvFile)
	if (inputCsvError is not None):
		return f"input CSV file error : {inputCsvError}"
	tumour = inputCsvResult[0]
	normal = inputCsvResult[1]

	tumour_tasks = seperateTask(tumour, parallelTask)
	normal_tasks = seperateTask(normal, parallelTask)

	for task in normal_tasks:
		inputPath = []
		outputFile = []
		samples = []
		for sample in task:
			inputPath.append(normal[sample])
			outputFile.append(f"{outputPath}/{sample}N_output.txt") 
			samples.append(sample)

		with ProcessPoolExecutor() as executor:
			tasks  = [executor.submit(featureCount, inputPath, outputFile, threadNo, outputPath) for inputPath, outputFile, samples in zip(inputPath, outputFile, samples)]

	for task in tumour_tasks:
		inputPath = []
		outputFile = []
		samples = []
		for sample in task:
			inputPath.append(tumour[sample])
			outputFile.append(f"{outputPath}/{sample}T_output.txt") 
			samples.append(sample)

		with ProcessPoolExecutor() as executor:
			tasks  = [executor.submit(featureCount, inputPath, outputFile, threadNo, outputPath) for inputPath, outputFile, samples in zip(inputPath, outputFile, samples)]

	resultDirCSV(tumour, 'T', outputPath)
	resultDirCSV(normal, 'N', outputPath)

	return None

def inputFileHandle(inputCsvFile):
	try : 
		input_csv = pd.read_csv(inputCsvFile , header=None)
		bam_list = input_csv.values.tolist()

		tumour = {}
		normal = {}

		for sample in bam_list:

			## check the format of csv
			if (len(sample) != 3):
				print("non valid input")
				break
			sample_id = sample[0]
			grouping = sample[1]
			route = sample[2]
			if (grouping == 'T'):
				tumour[sample_id] = route
			elif (grouping == 'N'):
				normal[sample_id] = route

		## check if normal and tumour sample are corresponding
		if (list(normal.keys()) != list(tumour.keys())):
			print("tumour and normal sample are not corresponding")

		return [tumour, normal], None

	except pd.errors.ParserError as err: 
		return None, err

def seperateTask(sampleDictionary, parallelTask):
	task_lists = []
	task_list = []
	n = 0

	for sample in sampleDictionary:
	    n += 1
	    if (n <= int(parallelTask)):
	        task_list.append(sample)
	    else:
	        task_lists.append(task_list)
	        task_list = []
	        n = 1
	        task_list.append(sample)

	task_lists.append(task_list)

	return task_lists

def featureCount(inputPath, outputFile, thread, outputPath):
	featureCountsCommand = f"featureCounts -p --countReadPairs -F GTF -t HERV -T {thread} -a {outputPath}/quantification.gtf -o {outputFile}  {inputPath}"
	try:
		subprocess.run(featureCountsCommand, shell = True, check = True)
	except subprocess.CalledProcessError as err:
		print(f"Error when running featureCounts : {err} \n")

def resultDirCSV(sampleDictionary, typeOfSample, outputPath):
	resultDir = []

	## add each sample to the dictionary
	for sample in sampleDictionary:
		resultDir.append({
			'sample name' : sample, 
			'type' : typeOfSample, 
			'result directory' : f"{outputPath}/{sample}{typeOfSample}_output.txt"
		})
		resultDir.append({
			'sample name' : sample, 
			'type' : typeOfSample, 
			'result directory' : f"{outputPath}/{sample}{typeOfSample}_output.txt.summary"
		})
		
	resultDir_dataFrame = pd.DataFrame(resultDir)
	resultDir_dataFrame.to_csv(f"{outputPath}/resultDirectories_{typeOfSample}.csv", index = False)







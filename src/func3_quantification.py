import subprocess
import pandas as pd
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

def func3_quantification(inputCsvFile, annotationFile, outputPath, threadNo, parallelTask, multiMapped):
	"""
	Use FeatureCounts to quantify each BLAST result remaining after filtering.
	Input : 
		inputCsvFile : path to the csv file with paths to the bam files
		annotationFile : the path of quantification.gtf
		outputPath : path for the output file
		threadNo : no. of thread to be used
		parallelTask : integer for the parallel task
	Output : resultDirectories_{T/N}.csv which record the path of the output of FeatureCounts
	"""
	# to check the process for server
	with open(f"{outputPath}/streamlit_quantification_process.txt", "w") as f:
		f.write("0")
	# get bam file directories from csv input
	inputCsvResult, inputCsvError = inputFileHandle(inputCsvFile)
	if (inputCsvError is not None):
		raise ValueError(inputCsvError)
		sys.exit(1)
	tumour = inputCsvResult[0]
	normal = inputCsvResult[1]

	# separate tasks into parallel list
	if len(tumour.keys()) != 0:
		tumour_tasks, tumourTasksError = separateTask(tumour, parallelTask)
		if (tumourTasksError is not None):
			raise ValueError(tumourTasksError)
			sys.exit(1)
		# Run parallel task
		paraTumourTaskError = processParallelTasks(tumour, tumour_tasks, "T", outputPath, threadNo, annotationFile, multiMapped)
		if (paraTumourTaskError is not None):
			raise ValueError(paraTumourTaskError)
			sys.exit(1)
		resultDirCSVerror_T = resultDirCSV(tumour, 'T', outputPath)
		if (resultDirCSVerror_T is not None):
			raise ValueError(resultDirCSVerror_T)
			sys.exit(1)

	if len(normal.keys()) != 0:
		normal_tasks, normalTasksError = separateTask(normal, parallelTask)
		if (normalTasksError is not None):
			raise ValueError(normalTasksError)
			sys.exit(1)
		paraNormalTaskError = processParallelTasks(normal, normal_tasks, "N", outputPath, threadNo, annotationFile, multiMapped)
		if (paraNormalTaskError is not None):
			raise ValueError(paraNormalTaskError)
			sys.exit(1)
		resultDirCSVerror_N = resultDirCSV(normal, 'N', outputPath)
		if (resultDirCSVerror_N is not None):
			raise ValueError(resultDirCSVerror_N)
			sys.exit(1)
		
	with open(f"{outputPath}/streamlit_quantification_process.txt", "w") as f:
		f.write("1")

	return None

def inputFileHandle(inputCsvFile):
	"""
	Handle the input csv file containing paths to the bam files of the samples
	Input : 
		inputCsvFile : csv file with paths to the bam files
	Output : list of tumour and normal dictionary with {sample_id : path}
	Return : output_list, error
	"""
	try : 
		input_csv = pd.read_csv(inputCsvFile , header=None)
		bam_list = input_csv.values.tolist()

		tumour = {}
		normal = {}
		samples = {}

		for sample in bam_list:

			# check the format of csv
			if (len(sample) != 3):
				raise ValueError("non valid csv input")

			sample_id = sample[0]
			grouping = sample[1]
			route = sample[2]
			if (grouping == 'T'):
				tumour[sample_id] = route
			elif (grouping == 'N'):
				normal[sample_id] = route

		# check if normal and tumour sample are corresponding
		# if (list(normal.keys()) != list(tumour.keys())):
		# 	raise ValueError("tumour and normal sample are not corresponding")

		return [tumour, normal], None

	except pd.errors.ParserError as err: 
		print(f"Error: please check yout csv input: {err}", file=sys.stderr)
		return None, err
	except ValueError as err:
		print(f"Error: please check your csv input: {err}", file=sys.stderr)
		return None, err 
	except Exception as err:
		print(f"Error: {err}", file=sys.stderr)
	
def separateTask(sampleDictionary, parallelTask):
	"""
	Separate task into parallel task list 
	Input : 
		sampleDictionary : list of dictionary for tumor and normal samples
		parallelTask : integer for the parallel task
	Return : 
		task_lists : list of list of task to run in parallel
	"""
	
	try:
		if (type(int(parallelTask)) != int):
			raise ValueError("number of parallel task has to be integer")
		
		task_lists = [[] for _ in range(int(parallelTask))]
		i = 0

		for sample in sampleDictionary:
			if (i < int(parallelTask)):
				task_lists[i].append(sample)
				i += 1
			else:
				i = 0
				task_lists[i].append(sample)
				i += 1

		# task_lists = []
		# task_list = []
		# n = 0

		# for sample in sampleDictionary:
		# 	n += 1
		# 	if (n <= int(parallelTask)):
		# 		task_list.append(sample)
		# 	else:
		# 		task_lists.append(task_list)
		# 		task_list = []
		# 		n = 1
		# 		task_list.append(sample)

		# task_lists.append(task_list)

		return task_lists, None
	except ValueError as err:
		print(f"Error: {err}", file=sys.stderr)
		return None, err 
	except Exception as err:
		print(f"Error: {err}", file=sys.stderr)
		return None, err 

def processParallelTasks(sampleDictionary, tasks, typeOfSample, outputPath, threadNo, annotationFile, multiMapped):
	"""
	Process parallel task using ProcessPoolExecutor
	Input : 
		sampleDictionary : the dictionary contain list of sampleID:samplePath
		tasks : list of task (list of sampleID)
		typeOfSample : T/N for tumour or normal sample
		outputPath : path for the output file
		threadNo : no. of thread to be used
		annotationFile : the path of quantification.gtf
	Output :
		results of featureCount
	"""
	try:
		with ProcessPoolExecutor() as executor:
			futures = []
			inputPaths = []
			outputFiles = []
			for task in tasks:
				for sample in task:
					inputPaths.append(sampleDictionary[sample])
					outputFiles.append(f"{outputPath}/{sample}{typeOfSample}_output.txt") 
				# Submit task to the executor
				futures.extend(executor.submit(featureCount, inputPath, outputFile, threadNo, annotationFile, multiMapped) for inputPath, outputFile in zip(inputPaths, outputFiles))

			# Handle the tasks
			for future in as_completed(futures):
				future.result() 
	except ValueError as err:
		print(f"Task error: {err}", file=sys.stderr)
		return err 
	except Exception as err:
		print(f"Task error: {err}", file=sys.stderr)
		return err 
	
	return None

def featureCount(inputPath, outputFile, threadNo, annotationFile, multiMapped):
	if multiMapped == "1":
		featureCountsCommand = f"featureCounts -M --fraction -p --countReadPairs -F GTF -t HERV -T {threadNo} -a {annotationFile} -o {outputFile} {inputPath}"
	else:
		featureCountsCommand = f"featureCounts -p --countReadPairs -F GTF -t HERV -T {threadNo} -a {annotationFile} -o {outputFile} {inputPath}"

	try:
		subprocess.run(featureCountsCommand, shell = True, check = True)
	except subprocess.CalledProcessError as err:
		raise ValueError(f"Error when running featureCounts: {err}")

def resultDirCSV(sampleDictionary, typeOfSample, outputPath):
	"""
	Produce the csv file containing the paths of the output files from featureCount
	Input : 
		sampleDictionary : the dictionary contain list of sampleID:samplePath
		typeOfSample : T/N for tumour or normal sample
		outputPath : path for the output file
	Output : resultDirectories_{T/N}.csv
	"""
	try:
		resultDir = []

		# add each sample to the dictionary
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
	
	except Exception as err:
		print(f"Error: {err}", file=sys.stderr)
		return err 
	
	return None




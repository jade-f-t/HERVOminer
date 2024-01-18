# HERVOminer

## Description
HERVOminer is a tool that capable of determining whether the query peptides originated from Human Endogenous retroviruses (HERVs) regions, quantify their expression and visualize the result.

## Installation

### System Requirements
HERVOminer is designed for Linux operating systems, full compatibility is only guaranteed on Linux.

### Installing Miniconda

HERVOminer utilizes tools from Bioconda, a channel for bioinformatics software of Conda. To ensure easy installation of the dependencies, it is recommended to install Miniconda if you do not have Conda or Miniconda installed. 

Please visit the [Miniconda installation page](https://docs.conda.io/projects/miniconda/en/latest/miniconda-install.html) and follow the installation instructions to install Miniconda.

### Setting Up a Conda Environment
It is recommended to create a Conda environment for using HERVOminer to avoid conflicts. 

Please use the following commands to create and activate a new environment:

```bash
conda create -n HERVOminer_env python=3.11
conda activate HERVOminer_env
```

### Installing HERVOminer

1. Clone the repository
```bash
git clone https://github.com/jade-f-t/HERVOminer.git
cd HERVOminer
```

2. Installing Dependencies \
Please use the following commands to install the dependencies:

```bash
conda install bioconda::subread=2.0.6
conda install bioconda::bedtools=2.31.0
pip install -r requirements.txt
```

3. Decompress Data \
Please decompress the included compressed data files before use, as they contain essential data required for the operation of HERVOminer.

```bash
cd data 
tar -xzvf HERV_ORF_dict.tar.gz
tar -xzvf HervOrfBlastpDB.tar.gz
cd ..
```

4. Grant execution permission to the HERVOminer.py script
```bash
chmod +x ./scripts/HERVOminer.py
```

## Implementation

HERVOminer is designed to determine whether the query peptides are derived from the HERV regions and quantify the corresponding sequences. HERVOminer processes 2 databases and 2 user inputs in 4 steps and generates 3 tables and 3 graphs of the analytical results. The 4 steps can be executed individually using 4 subcommands, or can be performed all at once with a single command.

### Input files format
1. RNA-seq CSV
- Description : \
A CSV file summarizing directories of post-processed alignment files for the samples.
Each tumor sample need to have its corresponding normal sample. The sample id fields of the tumor and corresponding normal sample have to be same and the sample type field should state whether the sample is tumor (T) or normal (N) sample. The directory has to be absolute directory.
- Format : \
3 fields : sample id, sample type (T / N), directory
- Example input format : 
| <!-- -->    | <!-- -->    | <!-- -->    |
|-------------|-------------|-------------|
| 10   | T   | /Users/jade-f-t/data/10T.bam   | 
| 10   | N   | /Users/jade-f-t/data/10N.bam   | 
| 21   | T   | /Users/jade-f-t/data/21T.bam   |
| 21   | N   | /Users/jade-f-t/data/21N.bam   |
| ...     | ...     | ...     |


2. query peptide 
- Description : \
The absolute directory of the query peptides FASTA file.
- Example input format : \
/Users/jade-f-t/data/input_peptide.fasta

### Subcommands

1. Protein BLAST Similarity Analysis (BLAST)
- Description : \
Find regions of similarity between query peptides and the open reading frames of HERV regions. 
- Usage :
```bash
./scripts/HERVOminer.py blastORF \
-p <path to the input peptide file> \
-o <path to the output file>
```
- Output : blastp_output.txt

2. BLAST results summarisation, Annotation file generation
- Description : \
Filter BLAST results of protein BLAST , and summarize their information.Then, generate the corresponding annotation file for the following quantification process.
- Usage :
```bash
./scripts/HERVOminer.py summarizeAnnotate \
-i <path to the input blastp output file get from Step 1> \
-o <path to the output file>
```
- Output : output_dict.json, quantification.gtf

3. Quantification (featureCounts)
- Description : \
Use FeatureCountsto quantify each BLAST result remaining after filtering.
- Usage :
```bash
./scripts/HERVOminer.py quantification \
-i <path to the input csv file with all of the bam file path> \
-a <path to the annotation file created by Step 2> \
-o <path to the output file> \
-t <no of threads to use for one sample> \
-n <number of sample to be analyzed in parallel computing>
```
- Output : output files of featureCounts, csv files for tumor and normal featureCounts   
output directories (resultDirectories_T.csv and resultDirectories_N.csv)

4. Data summarisation and visualization
- Description : \
Generate three tables and three groups of graphs from the quantification results.
- Usage :
```bash
./scripts/HERVOminer.py outputResult \
-p <path to the input peptide file> \ 
-i <path to the input csv file with all of the bam file path> \ 
-T <path to the resultDirectories_T.csv file> \ 
-N <path to the resultDirectories_N.csv file> \ 
-z <include region with zero count or not, 1 : with zero, 0: without zero> \ 
-o <path to the output file>
-d <set the dpi of the figures>
```
- Output : 3 Tables , 3 Figures





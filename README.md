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

2. Installing Dependencies
Please use the following commands to install the dependencies:

```bash
conda install bioconda::subread=2.0.6
conda install bioconda::bedtools=2.31.0
pip install -r requirements.txt
```

## Implementation
```bash
chmod +x ./scripts/HERVOminer.py
```
HERVOminer is separated into 4 subcommands : 





# Genomic Sequence Analysis Pipeline

This project provides a pipeline for genomic sequence analysis using Python and Bash scripts, utilizing Bowtie2 for sequence alignment, SAM/BAM file processing, and gene enrichment analysis.

## Table of Contents
- Overview
- Repository Structure
- Prerequisites
- Installation
- Usage
  - Running the Pipeline
  - Output
- Example Data
- License
- Contact

## Overview
The pipeline is designed to:
- Extract sequences from FASTQ files based on specific repeat patterns.
- Perform sequence alignment using Bowtie2.
- Analyze gene-level enrichment and differential expression between experimental and control samples.

The `dual-gRNA.py` script handles the main data processing, while `dual-gRNA.sh` is a Bash script to manage the overall workflow on an HPC cluster using SLURM and GNU Parallel.

## Repository Structure
```
/project_root
├── README.md               # Project documentation
├── .gitignore              # Files to be ignored by Git
├── requirements.txt        # Python dependencies
├── environment.yml         # Conda environment configuration
├── dual-gRNA.py             # Main Python script for analysis
├── dual-gRNA.sh             # Bash script for running the pipeline
├── data/                   # Example input data (provide your own)
└── results/                # Output results
```

## Prerequisites
- Python 3.9.18
- Conda (for managing dependencies)
- SLURM (for job scheduling on HPC clusters)
- Bowtie2 (for sequence alignment）
- GNU Parallel (for executing jobs in parallel)

## Installation

### Using Conda
1. Install Conda if not already installed. Follow the instructions at [Miniconda installation](https://docs.conda.io/en/latest/miniconda.html).
2. Create and activate the Conda environment:
   ```bash
   conda env create -f environment.yml
   conda activate sun
   ```

### Using pip
Alternatively, you can install the dependencies using `pip`:
```bash
pip install -r requirements.txt
```

## Usage

### Running the Pipeline
1. Modify the Bash script: Open `dual-gRNA.sh` and set the parameters, such as:
   - `samples_control` and `samples_experiment`: These arrays define the samples for the control group and the experimental group. By default, each array contains two samples (`"Control1"`, `"Control2"` for the control group and `"Experiment1"`, `"Experiment2"` for the experimental group). You can add more samples to these arrays as needed. For example, if you have 4 control samples and 4 experimental samples, you can modify them as follows:
     ```bash
     samples_control=("Control1" "Control2" "Control3" "Control4")
     samples_experiment=("Experiment1" "Experiment2" "Experiment3" "Experiment4")
     ```
     The script will process these samples in parallel.
   - `repeat` for the repeat sequence.
   - `genome_fa` for the reference genome file.

2. Submit the job using SLURM:
   ```bash
   sbatch dual-gRNA.sh
   ```

The script will:
- Extract sequences based on the repeat pattern.
- Perform alignment using Bowtie2.
- Process the SAM file and analyze gene enrichment for each sample.

### Output
- The analysis results will be saved in `/hpcfs/fhome/sunlet/dual`, as specified in the `dual-gRNA.sh` script.
- You will receive job completion or error notifications via email (configured in the Bash script).

## Example Data
- Place example FASTQ files in the `data/` directory to test the pipeline.

# License
This project is licensed under the Creative Commons Attribution-NonCommercial 4.0 International License.

## Contact
For questions or issues, please contact Sunflower_l@outlook.com.

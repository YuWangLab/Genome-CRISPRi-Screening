#!/bin/bash
#SBATCH --output=/hpcfs/fhome/sunlet/single/job.%j.out  # Directs the output of the job to a specific file, using job ID
#SBATCH --error=/hpcfs/fhome/sunlet/single/job.%j.err   # Directs the error output of the job to a specific file, using job ID
#SBATCH -N 1 -n 16                                                     # Requests 1 node and 16 cores
#SBATCH -J test                                                        # Names the job as "test"
#SBATCH --mem=100G                                                     # Allocates 100 GB of memory for the job
#SBATCH --partition=qcpu_23i                                           # Specifies the partition (queue) the job will run on
#SBATCH --mail-type=END,FAIL                                           # Specifies to send email notifications on job END and FAIL
#SBATCH --mail-user=Sunflower_l@outlook.com                            # Specifies the email address for sending notifications

mkdir -p /hpcfs/fhome/sunlet/single                                    # Creates the directory if it doesn't exist
cd /hpcfs/fhome/sunlet/single                                          # Changes the working directory to the specified path
source /etc/profile                                                    # Sources the global profile script for environment settings
conda activate sun                                                     # Activates the Conda environment named "sun"

# Define two arrays of samples, one for control and one for the experiment
samples_control=("Control1" "Control2")                                
samples_experiment=("Experiment1" "Experiment2")                       
repeat="ATCTACAACAGTAGAAATTC"                                          # Define a repeat sequence to be used in the script
genome_fa="library.fasta"                                              # Define the genome fasta file name

samples=("${samples_control[@]}" "${samples_experiment[@]}")           # Combine both sample arrays into one

# Use GNU Parallel to run a Python script in parallel across the defined samples
parallel 'python 1.fq_sam.py \
  --input_fastq1 {1}_1.fq \
  --input_fastq2 {1}_2.fq \
  --genome_fa '"${genome_fa}"' \
  --repeat '"${repeat}"' \
  --samples '"${samples[*]}"' \
  --samples_control '"${samples_control[*]}"' \
  --samples_experiment '"${samples_experiment[*]}"'' ::: "${samples[@]}"

echo "All samples processed."                                          # Echo a message indicating all samples have been processednull

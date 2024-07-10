#!/bin/bash

#SBATCH --job-name=MCME_Even_Butter_Climate
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1

module load miniconda

conda activate /home/wp288/.conda/envs/test_env
python ./experiments-folder/MCME-Even-Butter.py

#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --job-name=seq_races
#SBATCH --nodes=1
###SBATCH --ntasks=1
#SBATCH --array=1-30
#SBARCH --cpus-per-task=6
#SBATCH --mem 30g
#SBATCH --time=24:00:00
#SBATCH --output=seq_%a.out

sim=${SLURM_ARRAY_TASK_ID}
#sim=$(($s + 10))
echo $sim

module load R
Rscript 3_sim_seq.R --simulation "sim_${sim}"
#Rscript 4_smoothing.R --simulation "sim_${sim}"

#!/bin/bash
#SBATCH --partition=GENOA
#SBATCH --job-name=seq_ProCESS_locate
#SBATCH --nodes=1
#SBATCH --array=1-30
####SBARCH --tasks-per-node=1
#SBARCH --cpus-per-task=2
#SBATCH --mem 50g
#SBATCH --time=48:00:00
#SBATCH --output=seq_%a.out
####SBATCH --output=sim.out

sim=${SLURM_ARRAY_TASK_ID}
echo $sim

module load singularity
image="/orfeo/cephfs/scratch/cdslab/shared/containers/singularity/sarek_tumourevo/lvaleriani-process_validation-v1.img"
base="/orfeo/scratch/area/lvaleriani/utils_locate/simulations_rRACES"

singularity exec --bind /orfeo:/orfeo --no-home $image Rscript $base/4_plot.R --simulation "sim_${sim}"

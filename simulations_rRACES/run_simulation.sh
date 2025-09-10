#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --job-name=sim_ProCESS_locate
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBARCH --cpus-per-task=2
#SBATCH --mem 100g
#SBATCH --time=48:00:00
#SBATCH --output=sim.out


module load singularity
image="/orfeo/cephfs/scratch/cdslab/shared/containers/singularity/sarek_tumourevo/lvaleriani-process_validation-v1.img"
base="/orfeo/scratch/area/lvaleriani/utils_locate/simulations_rRACES"

#singularity exec --bind /orfeo:/orfeo --no-home $image Rscript $base/1_sim_tissue.R
singularity exec --bind /orfeo:/orfeo --no-home $image Rscript $base/2_put_mutations.R

#!/bin/sh

#SBATCH --job-name=ONT_QC
#SBATCH --partition=THIN
#SBATCH --nodes=1
###SBATCH --array=1-4
#SBATCH --tasks-per-node=1 
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=4:00:00
#SBATCH -A lade
#SBATCH --output=qc_%a.out

export PATH="/orfeo/scratch/area/lvaleriani/myconda/bin:$PATH"
source activate lr

sample_sheet='/orfeo/LTS/LADE/LT_storage/lvaleriani/scolorina_analysis/samplesheet.tsv'

idx=7
####$(($SLURM_ARRAY_TASK_ID+2))
sample=$(awk "NR==$idx{print \$1; exit}" $sample_sheet)
path=$(awk "NR==$idx{print \$2; exit}" $sample_sheet)
summary=$(find "$path" -type f -name "sequencing_summary_*" | head -n 1)

align="/orfeo/LTS/LADE/LT_storage/lvaleriani/scolorina_analysis/results/"
outdir="/orfeo/LTS/LADE/LT_storage/lvaleriani/scolorina_analysis/results/pycoQC"

pycoQC -f ${summary} -a ${align}/${sample}/fastcat.sorted.aligned.bam -o ${outdir}/${sample}_pycoQC_output.html

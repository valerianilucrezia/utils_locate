#!/bin/sh
#SBATCH --job-name=ONT_fastcat
#SBATCH --partition=THIN
#SBATCH --nodes=1
#SBATCH --array=1-6
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G
#SBATCH --time=5:00:00
#SBATCH -A lade
#SBATCH --output=fastcat_%a.out

export PATH="/orfeo/scratch/area/lvaleriani/myconda/bin:$PATH"
source activate fastcat

sample_sheet='/orfeo/LTS/LADE/LT_storage/lvaleriani/scolorina_analysis/samplesheet.tsv'

idx=$(($SLURM_ARRAY_TASK_ID+1))
sample=$(awk "NR==$idx{print \$1; exit}" $sample_sheet)
path=$(awk "NR==$idx{print \$2; exit}" $sample_sheet)

out=${path}/fastq_pass/fastcat
mkdir -p $out
rm -r $out/histogram/

echo $sample

fastcat \
    -s $sample \
    -f ${out}/per-file-stats.tsv \
    -i ${out}/per-file-runids.txt \
    --histograms ${out}/histogram/ \
    ${path}/fastq_pass/ > ${out}/${sample}.fastq

bgzip -@ 24 -f ${out}/${sample}.fastq 

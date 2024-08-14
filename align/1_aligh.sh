#!/bin/sh

#SBATCH --job-name=ONT_align
#SBATCH --partition=THIN
#SBATCH --nodes=1
####SBATCH --array=1-5
#SBATCH --tasks-per-node=1 
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=72:00:00
#SBATCH -A lade
#SBATCH --output=align_%a.out


module load singularity/3.10.4
module load java

sample_sheet='/orfeo/LTS/LADE/LT_storage/lvaleriani/scolorina_analysis/samplesheet.tsv'
#idx=$(($SLURM_ARRAY_TASK_ID+2))

for idx in 3 4 5 6 7
do

	sample=$(awk "NR==$idx{print \$1; exit}" $sample_sheet)
	path=$(awk "NR==$idx{print \$2; exit}" $sample_sheet)

	ref='/orfeo/LTS/LADE/LT_storage/lvaleriani/CNA/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'
	outdir="/orfeo/LTS/LADE/LT_storage/lvaleriani/scolorina_analysis/results/${sample}/"
	mkdir -p $outdir

/orfeo/LTS/LADE/LT_storage/lvaleriani/nextflow/nextflow run wf-alignment/main.nf\
	--fastq ${path}/fastq_pass/fastcat/ \
	--references $ref \
	--out_dir $outdir \
	-profile singularity \
	--threads 48 \
 	-c align.config \
	--per_read_stats false \
	--depth_coverage true \
	--disable_ping true 
done

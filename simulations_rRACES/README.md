## Simulations with rRACES

Required packages:
- `rRACES`
- `dplyr`
- `tidyr`
- `ggplot2`
- `patchwork`
- `optparse`
- `cli`

How to run scripts:

`Rscript 1_sim_tissue.R --outdir [out_dir_path] --type clonal`

`Rscript 2_put_mutations.R --outdir [out_dir_path] --type clonal`

Inside `2_put_mutations.R` there are some hardcoded parameters for simulation cobinations that can be changed.

For `3_sim_seq.R`, launch `run_sequencing.sh` where `#SBATCH --array=1-30` is the number of simulations you want to sequence.
The sequencing will be performed for coverage = [50,70,100] and purity = [0.3,0.6,0.9].



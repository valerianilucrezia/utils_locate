## Task
Run [HATCHET](#hatchet) and [Battenberg](#battenberg) copy number caller on:
- Organoids Nanopore data (use normal from Illumina)
- Organoids Illumina data
  

## Data location:
- Nanopore: `/orfeo/LTS/CDSLab/LT_storage/ORID0077_Scolorina`
```
└── PDO74
    └── 20240717_1531_2F_PAW71788_72fc2258
        └── run_5mCG_5hmCG
            ├── PDO74_5mCG_5hmCG_sorted.bam
            ├── PDO74_5mCG_5hmCG_sorted.bam.csi
└── PDO57_III
     └── 20240723_1454_2B_PAW71973_c9304556
        └── run_5mCG_5hmCG
            ├── PDO57_III_5mCG_5hmCG_sorted.bam
            ├── PDO57_III_5mCG_5hmCG_sorted.bam.csi
...
```

- Illumina: `/orfeo/LTS/LADE/LT_storage/lvaleriani/scolorina_analysis/short_read_data/bam`
>[!Caution]
>- `<sample_name>_A` normal sample
>- `<sample_name>` tumor sample
```
├── PDO55
|    ├── chr1-PDO55.bam
|    ├── chr1-PDO55.bam.bai
...
├── PDO55_A
├── PDO57_III
├── PDO57_III_A
├── PDO57_VII
├── PDO57_VII_A
├── PDO6
├── PDO61
├── PDO61_A
├── PDO6_A
├── PDO74
└── PDO74_A
```

## HATCHET
**Publication**: Zaccaria, Simone, and Benjamin J. Raphael. "Accurate quantification of copy-number aberrations and whole-genome duplications in multi-sample tumor sequencing data." Nature communications 11.1 (2020): 4301.

**Tool**: https://raphael-group.github.io/hatchet/#

## Battenberg
**Publication**: Nik-Zainal, S., Van Loo, P., Wedge, D. C., Alexandrov, L. B., Greenman, C. D., Lau, K. W., ... & Campbell, P. J. (2012). The life history of 21 breast cancers. Cell, 149(5), 994-1007.

**Tool**:
- https://github.com/Wedge-lab/battenberg
- https://github.com/cancerit/cgpBattenberg

## Output
- A folder for each sample (Illumina and Nanopore) with the output of both HATCHET and Battenberg.
- Scripts for the reprodubility of the results.

## Reference Files
- Nanopore: `/orfeo/LTS/LADE/LT_storage/lvaleriani/CNA/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna`
- Illumina: `/orfeo/LTS/CDSLab/LT_storage/ref_genomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta`



# Cancer Project Repository

## Contents
- [The structure of the project with a tree](#the-structure-of-the-project-with-a-tree)
- [Preprocessing of the RNAseq data](#preprocessing-of-the-rnaseq-data)


## The structure of the project with a tree

```
cancer_project/
├── data/
│   ├── databases/
│   ├── processed/
│   └── raw/
|
├── src/                            # source code for computational analysis       
│   ├── exploratory/                # exploratory analysis scripts
|
├── util/                           # utility scripts and helper functions for src 
|
├── ngs_preprocess/                 # RNAseq preprocessing pipeline
│   ├── 0_raw/                      # soft-linked FASTQ files from sra/GSE/parallel-fastq  
|   |
│   ├── 1_fastqc/                    # FastQC and MultiQC results pre- and post-trimming
│   │   ├── after_trimming_qc/
│   │   │   └── multiqc_data/
│   │   └── pre_trimming_qc/
│   │       └── multiqc_data/
|   |
│   ├── 2_trimmed/                  # trimmed FASTQ files
│   ├── 3_mapping/                  # STAR mapping results (BAM/SAM files)
│   ├── 4_qualimap/                 # Qualimap, RSeQC and other QC results of aligned reads  
│   ├── 5a_salmonout/               # Salmon quantification results from all data
│   ├── 5b_starout/                 # STAR-RSEM/HTSeq quantification results from all data
│   └── 6_final/                    # combined expression matrices from all data
|   |
│   ├── reference/                  # reference genomes, transcriptomes annotation files and indices 
│   │   ├── hg38_salmon/
│   │   └── hg38_star/
|   |
│   ├── scripts/                    # bash scripts used in the pipeline  
│   │   └── exploratory/            # exploratory scripts for testing 
|   |
│   └── sra/                        # directory where GEO are fetched/dumped
│       ├── GSE199451/              # specific GEO study
|           ├── parallel_fetch/     # output of parallel-fetch SRA files
|           ├── parallel_fastq/     # output of parallel-fastq SRA to FASTQ conversion
|           ├── SRR_Acc_List.       # list of SRA run accessions
|
|
|
└── results/
    ├── jpg/
    └── pdf/
```

### Preprocessing of the RNAseq data
The RNAseq data preprocessing pipeline is implemented in the `ngs_preprocess/` directory. The [README.md](ngs_preprocess/README.md) file inside that directory contains: <br> 
- An .md/html explaining each steps and how to create the directories.
- Detailed instructions on how to run the pipeline with a table explaining the scripts.
- Dependencies yml files for the miniconda environment.
- Descriptions of the reference genomes/transcriptomes and annotation files used.
- Comparisons of different tools during development e.g. alignment rates, quantification methods, running times, RAM/CPU usage, and QC metrics of some "bad" samples with low mapping percentages.
- The aim of the comparisons is to decide which tools and parameters to use/consider in a snakemake workflow. 

![image: pipeline](./results/jpg/pipeline.jpg)

<br>

**Could a Snakemake workflow be usefull here?** <br>
Running the pipeline above locally with limited CPU/RAM resourses does not warrant the use of a Snakemake workflow. The reason I made the comparisons of cores and RAM usage, is to decide whether it would be beneficial and feasible to run 2 samples in parallel with 11cores/32GB RAM locally. However, even *light tools* can be quite heavy on occasions. For instance, **running Salmon** with i) decoy, ii) ValidateMappings, iii) bootsrapping and iv) bias has **high memory footprint (>50GB RAM)** making RAM usage a limiting factor. Even if a workstation of 12-16 threads is good, depending on the analysis, a single sample could **skyrocket RAM usage**. Furthermore **FPPE samples have a higher rate of spurious reads** due to low quality of RNA making the abover options necessary. With 62GB i run successfully 2 heavy jobs for Salmon (all options enabled) in parallel, whereas in 3 some samples crashed. When running STAR, 2 jobs require far more RAM usage than what is available locally. <br>

Obviously in HPC clusters with 128GB+ RAM and 64+ cores, Snakemake would be a better option and running multiple samples with less cores would be more eficcient than running 1 sample as fast as possible. Furthermore, sceduling the heavy work during parallelization cpuld also be a wise option. Furthermore, the information about CPU/RAM usage and running times would be useful for anyone trying to run the pipeline in an HPC cluster with job schedulers like SLURM/SGE or Snakemake.

**Verdict based on local runs:** <br>
t's better to run two samples in parallel with half of your threads each rather than maxing out 1 sample with all threads. RAM usage should be monitored based on requirements of the analsysis and sample quality. <br>
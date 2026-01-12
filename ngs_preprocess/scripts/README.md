## All the bash scripts used for pre-processing the RNA-seq data are available in this directory.
- In the exploratory directory, there are some additional scripts that were used for testing and exploring different tools and methods.

In the box below is a list of the main scripts used in the pipeline along with a brief description of their functionality. The scripts are hyperlinked. 

| Script Name                | Description                                                  |
|----------------------------|--------------------------------------------------------------|
| [`process_fastqc`](./process_fastqc) | Performs quality control using FastQC and MultiQC on FASTQ files. The output is stored in `1_fastq` and uses GNU for parallelization. | 
[`trimm_fastp`](./trimm_fastp)   | Finds the soft-linked FASTQ files from `0_raw` and perform trimming using FastP. The output is trimmed FASTQ files in ``2_trimmed`` directory along with the HTML report. Uses GNU for parallelization. With **5 cores**, it takes around **1 min to filter ~60million** reads. The reads were filtered for polyG tails and trimmed based on quality in a sliding window of 4 starting from the `tail/3'end` with minimum score of phred=20. The script is set to **4 jobs** (20cores), and it takes ~80seconds for 4 pair-end samples. Then **FASTQC paralleled** is run for the new samples, followed by a multiQC. The FASTQC files are saved with a `_trimmed_1_fastq.gz` name tag in the `1_fastq` directory.|
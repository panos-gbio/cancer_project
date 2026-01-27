# Bash Notes for NGS Preprocessing

## Table of Contents

- [Setting up the Project and Reference files](#setting-up-the-project-and-reference-files)
- [Access FASTQ files with SRA-toolkit](#access-fastq-files-with-sra-toolkit)
    - [Installing SRA-toolkit](#installing-sra-toolkit)
    - [Download one SRR ID from GSE183947 using prefetch and fasterq-dump](#download-one-srr-id-from-gse183947-using-prefetch-and-fasterq-dump)
    - [Download multiple SRR IDs from GSE183947 in a text file.](#download-multiple-srr-ids-from-gse183947-in-a-text-file)
    - [Downloading FASTQ files without SRA-toolkit.](#downloading-fastq-files-without-sra-toolkit)
- [Make STAR/Salmon index for RNA-seq alignment](#make-starsalmon-index-for-rna-seq-alignment)
- [Quality control of the FASTQ files (FASTQC and MULTIQC)](#quality-control-of-the-fastq-files-fastqc-and-multiqc)
- [Trimming with FastP](#trimming-with-fastp)
- [Pseudo-alignment and Quantification with Salmon](#pseudo-alignment-and-quantification-with-salmon)
- [Mapping with STAR aligner](#mapping-to-genome-with-star-aligner)
- [Post Alignment QC (RSeQC, Qualimap, MultiQC)](#post-alignment-qc-rseqc-qualimap-multiqc)
- [Quantification with RSEM](#quantification-with-rsem)
- [Some technical details in GNU parallel](#some-technical-details-in-gnu-parallel)


## Setting up the Project and Reference files

First I will create the project folder structure in the cancer_project directory. Normally I will create a GSE specific subdirectory in the sra folder to store the SRA/FASTQ url list file. I will use either script 1a or 1b to fetch my FASTQ files in GSE specific folders. The rest of the folders are self explanatory.

The genomes from different databases (ENSEMBL, NCBI, UCSC) may have different names for the same chromosomes. It is important to keep track of the version of the genome used for alignment and quantification.

- ESNSEMBL: GRCh38.p14 V115, names the chromosomes as 1,2,3,...X,Y,MT
- UCSC: hg38, names the chromosomes as chr1, chr2,...chrX, chrY, chrM
- NCBI T2T: CHM13v2.0, names the chromosomes as chr1, chr2,...chrX, chrY, chrM plus other contigs. Also the existence of other contigs may affect the mapping and quantification, or the use of older tools that might not consider their naming. 

```bash
mkdir -p ${HOME}/cancer_project/ngs_preprocess/{0_raw, 1_fastqc 2_trimmed, 3_mapping, 4_qualimap, 5a_salmonout, 5b_starout, reference, scripts, sra}; mkdir reference/{hg38_star, hg38_salmon, hg38_rsem}

# make a symbolic link of sra data to 0_raw after fetching the sra data.
ln -s ${HOME}/cancer_project/ngs_preprocess/sra/GSE199451/parallel_fetch ${HOME}/cancer_project/ngs_preprocess/0_raw/

```

STAR requires a genome fasta file and a gtf annotation file to make an index for alignment to the genome. Salmon require a transcriptome annotation file for quantitation/(pseudo)alignment to the transcriptome. RSEM also requires to generate a reference to run the calculation functions. I need to got ENSEMBL and find the latest version of the human genome (GRCh38) and download both files. After I decompress in the reference dir I will use STAR to make the index and save the output to hg38_star.

- ENSEMBL links from human genome [release-115](https://ftp.ensembl.org/pub/release-115)
- Fasta file: https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz 
- GTF file: 	https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
- Transcriptome file: https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

I have also found the T2T version of the genome in NCBI ftp server:

- The link is here [https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/]. The NCBI version GTF file contains HUGO gene nomenclature that might affect downstream analysis since HUGO-multiple ENSEMBLE IDs. 
- The T2T version is not annoatted and updated as much as the GRCH38, which is used by most databases, aligners and diagnostic labs. 

I have found other resources of the T2T with all the annotation files

- The github project [link](https://github.com/marbl/CHM13)
- The UCSC genome server of the T2T-CHM13 v2.0/hs1 and other annotations [link](https://hgdownload.soe.ucsc.edu/downloads.html#human)

```bash
# download references 
wget -P ${HOME}/cancer_project/ngs_preprocess/reference https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget -P ${HOME}/cancer_project/ngs_preprocess/reference https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
wget -P ${HOME}/cancer_project/ngs_preprocess/reference https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

# the we decompress the files 
gunzip ${HOME}/cancer_project/ngs_preprocess/reference/*.gz

# finally we can use find function to check where the reference is located in cancer_project dir
find  ${HOME}/cancer_project" \
  -type f \
  -name 'Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz' \
  -print
```
The arguments to make the STAR index are in the script below. The read length of the specific NGS library were 150bp pair end so the overhang distance will be 149. The argument is used to build the splice junction database (how much sequence around annotated junctions is used). <br>
**Important** In wsl2 not all computer's RAM is allocated so I need to check how much is available. For STAR, we need high RAM (>30GB) and enough CPU cores to speed the process up during indexing. 

```bash
# check RAM and CPU availability in WSL2-Ubuntu
free -h
nproc 
```
I may also need to change the available RAM allocated to WSL2 using powershell by creating in `~/.wslconfig`. The script should be written in powershell and saved to the windows `C:\Users\yourusername\.wslconfig` path. 

```powershell
[wsl2]
# Max RAM WSL2 VM can use (default is 50% of Windows RAM) :contentReference[oaicite:3]{index=3}
# Leave headroom for Windows (browser, IDE, etc.)
memory=48GB

# Max logical processors WSL2 can use (default is all logical processors) :contentReference[oaicite:4]{index=4}
processors=24

# Swap is disk-based memory used when RAM is exhausted (default 25% of memory) :contentReference[oaicite:5]{index=5}
# Setting some swap helps prevent OOM-kills, but it will be slower than RAM.
swap=10GB
```

Some aliqnment QC tools such as RSeQC require BED files for gene annotation. I can either convert the GTF to BED using ucsc functions or download pre-made BED files from RSeQC website [here](https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/)

Below I will convert the .gtf file to .bed file, which will be used in QC section. 


```bash

# install and check whereabouts 
conda install bioconda:ucsc-genepredtobed
conda install bioconda:ucsc-gtftogenepred
which gtfToGenePred
which genepredToBed

# convert in two steps
gtfToGenePred \
"${HOME}/cancer_project/ngs_preprocess/reference/Homo_sapiens.GRCh38.115.gtf" \
"${HOME}/cancer_project/ngs_preprocess/reference/Homo_sapiens.GRCh38.115.genePred"

genePredToBed \
"${HOME}/cancer_project/ngs_preprocess/reference/Homo_sapiens.GRCh38.115.genePred" \
"${HOME}/cancer_project/ngs_preprocess/reference/Homo_sapiens.GRCh38.115.bed"

# remove intermediate file
rm -r "${HOME}/cancer_project/ngs_preprocess/reference/Homo_sapiens.GRCh38.115.genePred"

```
A final step instead of running all the bed file for some slow functions of RSeQC, is to either generate the bed annotation of housekeeping gene of my current version or download the one that  is available [here](https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/).

```bash 
wget 

```


## Access FASTQ files with SRA-toolkit 

The SRA archives a vast amount of sequerncing data. There are some prefixes in a typical SRA accession number: '

- SRP#: Study
- SRX#: Experiment
- SRS#: Sample
- SRR#: Run

Furthermore a typical SRA would probably have a series matrix file (usually in .txt format) that contains the metadata of all the samples in the study. It could also contain the SRA accession number as well as the count matrix of the gene expression data/pre-processed data. Information on samples, sequencing instruments, adapters used, etc. can be found in the series matrix file. <br>
Other things available in SRA, could be some normalized data files () in GEO supplementary files.

Where can I find SRA data -> Gene expression omnibus (GEO) database: https://www.ncbi.nlm.nih.gov/geo/



### Installing SRA-toolkit 

I cant always use wget or curl to download SRA files directly form GEO. Maybe in other databases such as ENA where FASTQ files are directly accessible via ftp. For GEO, I need to use the SRA-toolkit to download the SRA files first, unzip, then convert them to FASTQ files. Very good instructions in this [video](https://www.youtube.com/watch?v=zE851fWCYQg). In the following example:

- I will use data from here [GSE183947](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183947).
- The SRA-toolkit can be downloaded from here: [https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit](https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit)
- Configuration instructions [here](https://github.com/ncbi/sra-tools/wiki/03.-Quick-Toolkit-Configuration)
- In the configuration set a tmp folder path for SRA. My path was `${HOME}/Downloads/tmp`.
- For prefetch and fasterq-dump instructions see [here](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump) (prefetch downloaded around 200MB per minute, a SRR of 3GB took around 15 minutes - these are ~30million reads RNAseq file).
- It is efficient to first use prefetch, create a directory with downloaded SRA files, then use fasterq-dump to convert them to FASTQ files.
- Obviously parallelization is a must for multiple SRR downloads/conversions. 


In SRA-run selector, I can select the runs I want to download. I can get the SRR IDs in to a text file for batch downloading using the SRA-toolkit. Of course once fetched, the data must be unzipped and splitted into Read1 and Read2 FASTQ files if it is paired-end sequencing data. Note that the **bin filder** of the SRA-toolkit contains the executables such as `fastq-dump` and `prefetch`. So either I use them by specifying the full path or I add the bin folder to the PATH variable in my `.bashrc` file. In the latter case, executables can be used directly from any folder in the terminal. 

SRA-toolkit for Ubuntu 22.04 LTS

```bash
# make a downloads directory and create tmp folder 
mkdir -p ~/Downloads/tmp; cd ~/Downloads

# download the tar file with wget in the Downloads folder
wget --content-disposition https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.2.1/sratoolkit.3.2.1-ubuntu64.tar.gz

# extract the tar file in the ~/Downloads folder
tar -xvzf sratoolkit.3.2.1-ubuntu64.tar.gz -C ~/Downloads

# configure using the executable in the bin folder, use i option 
./sratoolkit.3.2.1-ubuntu64/bin/vdb-config -i
```



### Download one SRR ID from GSE183947 using prefetch and fasterq-dump

The typical procedure is as follows:
- We can use the SRA by providing the full path to the bin folder or adding the bin folder to the PATH variable in the `.bashrc` file.
- `prefetch [SRR] -O [path]` downloads the SRA file to the specified path. For larger downloads, need to specify the file size limit using the `--max-size u` or `--max-size 100t` options. 
- for a txt of SRA IDs: `prefetch --option-file [file.txt] -O [path]` is used. 
- `fasterq-dump [SRR] -O [path] --split-files` converts the SRA file to FASTQ files and splits them if paired-end. The output FASTQ files are stored in the specified path.
- The commands work at CWD or at specified paths. Prefetch generates a folder based on the name of the SRR ID of that sample. Fasterq-dump needs to be run in the same directory where the SRA file was downloaded or it can be provided with the full path to the SRA file. Still i need an output directory to be specified if I dont want the FASTQ files to be created in the CWD.
- **SOS** a smart move could be to make a symbolic link to the directory where we want to save the FASTQ files. Use `ln -s [target file/folder] [link name]` to create a symbolic link, and use it as anrgument for -O option in fasterq-dump. 

```bash 
# I will make a folder for downloaded data in another directory called practice
mkdir -p ~/practice/sra_data

# first use prefetch to get the SRA files 
~/Downloads/sratoolkit.3.2.1-ubuntu64/bin/prefetch -O ~/practice/sra_data SRR15852393

# after prefetch, use fasterq-dump to convert SRA to FASTQ files in the sra_data
# it must be run in the same directory where the accession was downloaded
# here it is in the sra_data folder, but i will access it from practice
cd ~/practice
~/Downloads/sratoolkit.3.2.1-ubuntu64/bin/fasterq-dump --split-files ./sra_data/SRR15852393 -O ./sra_data

# we can also work with gzip compression to save space 
~/Downloads/sratoolkit.3.2.1-ubuntu64/bin/fasterq-dump --split-files ./sra_data/SRR15852393 -O ./sra_data --gzip

```



### Download multiple SRR IDs from GSE183947 in a text file.

First in GEO I will extract a text file containing 4 SRR IDs from the SRA run selector. The generated `SRR_Acc_List.txt` is saved at the `~/practice/sra_data` folder. Then I will use prefetch and fasterq-dump to download and convert them to FASTQ files.
- the prefetch will work like this `prefetch --option-file [file.txt] -O [path]`
- It took almost 42 minutes to prefetch 4 SRA files. 

```bash
~/Downloads/sratoolkit.3.2.1-ubuntu64/bin/prefetch -O ~/practice/sra_data --option-file ~/practice/sra_data/SRR_Acc_List.txt --max-size u

```

A more parrallelized approach since each sample could take around 10-30 minutes depending on the size. 

```bash
#!/usr/bin/env bash

# I will use GNU parallel to fetch multiple SRAs at once. 

JOBS=10 
SRA_FILE=~/practice/sra_data/SRR_Acc_List.txt
OUTPUT_FETCH=~/practice/sra_data/parallel_fetch
OUTPUT_FASTQ=~/practice/sra_data/parallel_fastq

# functions from SRA toolkit - paths to variables 
PREFETCH=~/Downloads/sratoolkit.3.2.1-ubuntu64/bin/prefetch
FASTERQ_DUMP=~/Downloads/sratoolkit.3.2.1-ubuntu64/bin/fasterq-dump

if [[ ! -s "${SRA_FILE}" ]]; then
    echo "SRA text file is empty or not found."
    exit 1
else
    echo "SRA text file list exists and there are $(cat "${SRA_FILE}" | wc -l) SRAs to fetch."
fi 

# check out the content of the SRA list file

while IFS= read -r read; do
    echo "Fetching ${read}..."
done < "${SRA_FILE}" # redirecting variable expansion 

# will use GNU parallel to fetch the SRAs, convert into a cat file and pipe
# cat retains line breaks and GNU recognises each line as a seperate input job
# I will use the available path 

cat "${SRA_FILE}" \
    | parallel -j "${JOBS}" --linebuffer --tag \
        "${PREFETCH}" -O "${OUTPUT_DIR}" --max-size u {}
```
The next step is to use fasterq-dump to convert the downloaded SRA files to FASTQ files. Again I will use GNU parallel to speed things up instead of converting the files one by one.  
- The fasterq-dump must be run in the directory where the SRA files were downloaded or provided with the full path to the SRA files. Also it uses by default 6 threads per job, so I will limit the number of parallel jobs to 3 at a time to avoid overloading the system.
- The output FASTQ files will be stored in a different folder.
- I can zip each FASTQ file using gzip commad as `gzip ./*.fastq` after the fasterq-dump is finished. Every FASTQ file will be zipped separately.

```bash
#!/usr/bin/env bash

# Continue with fasterq-dump parallelization here 
if [[ -z "$(ls -A "${OUTPUT_FASTQ}")" ]]; then

    echo "FASTQ directory is empty, starting fasterq-dump..."
    echo "Changing to fetch directory"
    cd "${OUTPUT_FETCH}" || { echo "Could not change directory"; exit 1; }

    # it outputs the filename without the preceding path using %P 
    find . -maxdepth 1 -type d -name "SRR*" -printf '%P\n' \
        | parallel -j 3 --linebuffer \
            --tagstring "Job {#} of ({})" \
            "${FASTERQ_DUMP}" -O "${OUTPUT_FASTQ}" --split-files -e 6 {}

else
    echo "FASTQ output directory is not empty, no need to run fasterq-dump. Remove anything non FASTQ."
fi

```

After some bad experience with `gzip` speed, a parallel implementation of gzip is also possible either with pre-existing tools such as `pigz` or using GNU parallel to gzip each FASTQ file concurrently. I have read that `pigz` is much faster than `gzip` since it uses multiple cores for compression.

```bash
# download and install pigz if not already installed
sudo apt-get install pigz

# use pigz for faster compression of FASTQ files to .gz
clear; echo "Compressing FASTQ files with pigz..."
pigz -p 12 "${OUTPUT_FASTQ}"/*.fastq

```

All the steps discussed above are saved in the same scrip at `${HOME}/cancer_project/ngs_preprocess/scripts/sra_parallel_fetch_fastq.sh` <br>


### Downloading FASTQ files without SRA-toolkit.

If the data is public in GEO, it may be accessible via ftp or http links in **ENA - European Nucleotide Archive**. For example, in ENA I can search for the same GSE183947 study using the ID and get direct links to the FASTQ files. Then I can use wget or curl to download them directly without using SRA-toolkit. This is a much simpler approach and I can parallelize the downloads using GNU-parallel which will speed things up. I am parallelizing the curl or wget proccesses using GNU-parallel as shown below.<br>

1) First I can use the SRA-explorer tool [here](https://sra-explorer.info/) that can also generate a bash script to download multiple FASTQ files and check whether they are available in EBI ftp or ENA servers.


2) Another option if prefetch is troublesome is to use wget to download the .sra file directly from the http link in the SRA-run selector. Then use fasterq-dump to convert it to FASTQ files. The SRA server is limited to 10 jobs for prefetch, maybe wget can bypass this. 

3) There is also a tool called fastq-dl that can download FASTQ files directly from ENA using the SRR IDs. The tool is available [here](https://github.com/rpetit3/fastq-dl)

**OTHER IDEAS - PARALLELIZATION in Fetching** <br>
Some Ideas in parallelizations and downloading multiple files or running multiple prefetch commands. Interesting discussion [here](https://github.com/ncbi/sra-tools/issues/560) <br>
GNU-parallel is also discussed [here](https://www.biostars.org/p/63816/) <br>
The SRA-explorer tool [here](https://sra-explorer.info/) can also generate a bash script to download multiple SRA files in parallel using prefetch and fasterq-dump.

### Make STAR/Salmon index for RNA-seq alignment 

The scripts below will be used for making the STAR and Salmon indices at `${HOME}/cancer_project/ngs_preprocess/reference` directory in the respective subfolders.

For STAR, the choice below will result in an index where STAR can map to transcriptome and genome simultaneously, and will select the best alignment. The gtf file is used to build the splice junction database and the option `--sjdbOverhang` will not be used during mapping but only during the genomeGenerate step. In the manual, genomeDir contents is described as genome sequence + suffix arrays + splice junction coordinates + transcript/gene information.

```bash
STAR --runThreadN 22 \
--runMode genomeGenerate \
--genomeDir ~/cancer_project/ngs_preprocess/reference/hg38_star \
--genomeFastaFiles ~/cancer_project/ngs_preprocess/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile ~/cancer_project/ngs_preprocess/reference/Homo_sapiens.GRCh38.115.gtf \
--sjdbOverhang 100 \

```
Salmon requires a transcriptome reference and perform quasi-mapping. Some important parameters for the creation of the Salmon index are the CPU cores and whether we want to generate a decoy transcriptome. There are some resources discussing why we need a decoy transcriptome to locate spurious reads that affect quantification:
- The official github page of salmon is [here](https://github.com/COMBINE-lab/salmon)
- A tutorial for salmon index and decoy transcriptome creation is [here](https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/)
- The salmon [paper](https://www.nature.com/articles/nmeth.4197) and the effect of spurious reads [here](https://www.biorxiv.org/content/10.1101/2021.01.17.426996v1)

**Important Note**
When I tried to install salmon in my existing conda env, the channel installed an old version of it and created dependency problems. It could be better to create a separate conda env for salmon only. The salmon index script I used is here: <br>

```bash
# for the decoy, I need to extract the names from the genome Fasta and 
# remove `>` from the names with sed 

grep -E "^>"  "${HOME}"/cancer_project/ngs_preprocess/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa | cut -d " " -f 1 | sed "s#>##" > decoy.txt 

cd "${HOME}/cancer_project/ngs_preprocess/reference"

# concatanate transcriptome and genome .fa > gentrome 
cat Homo_sapiens.GRCh38.cdna.all.fa Homo_sapiens.GRCh38.dna.primary_assembly.fa > gentrome.fa

# check whether decoy names are the sames as chromsome names in gentrome.fa
grep -E "^>1" gentrome.fa 

# before using Salmon make a salmon enviroment 
conda create -n salmon_env -c conda-forge -c bioconda "salmon=1.10.3"

# I encountered problems with salmon and my rest of the enviroment so I made conda env specifically
conda activate salmon_env
conda list | grep -i salmon 


# make the index with 12 cores 
salmon index -t gentrome.fa -d decoy.txt -p 12 -i hg38_salmon -k 31
echo $; 
```

Finally, I will probably need to create a RSEM reference index with GTF annotation file. The command below is followed when we already have a STAR index, using the same GTF and FASTA files/versions as STAR. 

```bash 
# need to make the directory 
mkdir -p -- "${HOME}/cancer_project/ngs_preprocess/reference/hg38_rsem"

# use GTF and Fasta genome files 
rsem-prepare-reference --gtf "${HOME}/cancer_project/ngs_preprocess/reference/Homo_sapiens.GRCh38.115.gtf" \
"${HOME}/cancer_project/ngs_preprocess/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
"${HOME}/cancer_project/ngs_preprocess/reference/hg38_rsem/hg38_rsem"

```

### Quality control of the FASTQ files (FASTQC and MULTIQC)

I will use directly the FastQC command over all the linked FASTQ.gz file in the input directory. A specific note here. The `&` indicates to run asynchronously in the background and shell does not wait to finish. So, this means that many Fastq files will be passed to the program. That i swhy we need to put a maximum 8 files per batch.

```bash
#!/usr/bin/env bash

INPUT=$HOME/cancer_project/ngs_preprocess/0_raw
OUT=$HOME/cancer_project/ngs_preprocess/1_fastq/
MAX_FILES=8

for fastq in "${INPUT}/*.fastq.gz"; do
     fastqc \
     -o "${OUT}" \
     $fastq &
     jobs -pr | wc -l
done

wait
```
Some technicalities in FASTQ
- Overepresented sequences and contaminations can influence the per base sequence content 


### Trimming with FastP 

FastP supports parallelization and runs with 3-4 cores to speed up the trimming process. Here I have some notes on the FastP tool commands 
- There are 2 inputs for pair-end data: `-i for Read1 and -I for Read2`
- The output trimmed files are specified with `-o` and `-O` for Read1 and Read2 respectively. The name of the output files can be specified directly.
- The `--length_required` or `-l` option sets the minimum length of reads to be retained after trimming. Usually set to 36bp.
- For adapter trimming, in paired-end data the `--detect_adapter_for_pe` option enables automatic detection and removal of adapters from both reads. It uses overlapping regions between Read1 and Read2 to identify adapter sequences. There is also the option to provide custom adapter sequences using `--adapter_sequence` and `--adapter_sequence_r2` for Read1 and Read2 respectively. We need to find the Ilumina adapter sequences used in the specific library prep kit (Truqseq adpaters, etc..)
- There is also the `--correction` option which enables base correction for paired-end data. It uses overlapping regions between Read1 and Read2 to correct sequencing errors.
- Another usefull option to delete polyG and polyX tails is `--trim_poly_g` and `--trim_poly_x` respectively.
- The `--cut_right` and `--cut_left` options enable cutting bases from the right and left ends of reads based on quality scores. The specific number of bases to cut can be specified with `--cut_right_mean_quality` and `--cut_left_mean_quality`. These defaults are set to 20 for quality and 4 for window size. These options start scanning from the respective ends. Then when a window fails the wuality threshold, it cuts the read from that window and so on. The window size is controlled by `--cut_right_window_size` and `--cut_left_window_size`.
- Also sliding window trimming can be performed with `--cut_window_size` and `--cut_mean_quality`. Finally `--cut_front` and `--cut_tail` can be used to trim low-quality bases from the front and tail of reads. The number of bases to trim can be specified with `--cut_front_mean_quality` and `--cut_tail_mean_quality`.
- There is a difference between `--cut_right` and `--cut_front/--cut_tail` options. The `--cut_right` option scans from the 5' end towards the 3' end and cuts bases once the quality threshold is not met. The rest of the read is cut off even if the quality improves later in the read. The `--cut_tail` option starts at the 3' end and trims low-quality bases from the end of the read until a base with quality above the threshold is encountered. The rest of the read is retained. The former is more agressive, possibly for variant calling, while the latter is more conservative, retaining more of the read.

```bash
# a part of a function in a bash script
fastp \
    -i "${f1}" \
    -I "${f2}" \
    -o "${OUT}/${BASENAME}_trimmed_1.fastq.gz" \
    -O "${OUT}/${BASENAME}_trimmed_2.fastq.gz" \
    --thread ${FASTP_THREADS} \
    --qualified_quality_phred 30 \
    --length_required 36 \
    --detect_adapter_for_pe \
    --correction \
    --trim_poly_g \
    --trim_poly_x \
    --cut_tail \
    --cut_tail_mean_quality 20 \
    --cut_tail_window_size 4 \
    --html "${OUT}/${BASENAME}_fastp.html"
    
    echo "Completed trimming for ${BASENAME}"

# an example of running fastp directly in terminal with just 2 samples

fastp \
    -i ${HOME}/cancer_project/ngs_preprocess/0_raw/SRR15852393_1.fastq.gz \
    -I ${HOME}/cancer_project/ngs_preprocess/0_raw/SRR15852393_2.fastq.gz \
    -o "${HOME}/cancer_project/ngs_preprocess/2_trimmed/SRR15852393__trimmed_1.fastq.gz" \
    -O "${HOME}/cancer_project/ngs_preprocess/2_trimmed/SRR15852393__trimmed_2.fastq.gz" \
    --thread 6 \
    --qualified_quality_phred 30 \
    --length_required 36 \
    --detect_adapter_for_pe \
    --correction \
    --trim_poly_g \
    --trim_poly_x \
    --cut_tail \
    --cut_tail_mean_quality 20 \
    --cut_tail_window_size 4 \
    --html "${HOME}/cancer_project/ngs_preprocess/2_trimmed/SRR15852393_fastp.html"
```

#### Options used in Trimmomatic (for comparison)

Trimmomatic is another tool used for trimming FASTQ files. This tool is slower than FastP and requires manual input of the adapter sequences. FastP offers the choice of manual input, however, it uses an algorith to locate the adapters in case there is no input. Some options used in Trimmomatic are:
The reading and trailing options remove bases at the start 5' end and at the 3' end of the reads before the index/adapter if teir quality is below a certain threshold. 
- `LEADING:3` : removes leading bases with quality below 3.
- `TRAILING:3` : removes trailing bases with quality below 3.
- `MINLEN:36` : drops reads that are shorter than 36 bases after trimming.
- `SLIDINGWINDOW:4:20` : performs sliding window trimming, cutting once the average quality within the window of 4bps falls below 20. The window size is 4 bases. For long reads this is useful to remove low-quality ends. So, many reads will be trimmed into smaller lengths and of those only the ones that are above 36bp will be retained.

### Pseudo-alignment and Quantification with Salmon

After trimming the FASTQ files, I will use Salmon to perform quasi-mapping and quantification of the reads to the transcriptome. The important options used in Salmon quant are:

A simple bash script for one paired-end sample is shown below. I also timed the process:

```bash
#!/usr/bin/env bash 

# set parameters for FASTQ and quant storing 
SALMON_INDEX="${HOME}"/cancer_project/ngs_preprocess/reference/hg38_salmon 
INPUT="${HOME}"/cancer_project/ngs_preprocess/2_trimmed
OUT="${HOME}"/cancer_project/ngs_preprocess/5a_salmonout
TIME_OUT="${HOME}"/cancer_project/ngs_preprocess/scripts/exploratory/time.log
SALMON_THREADS=7

SECONDS=0

# Test for 30 million paired end reads after trimming. 
# Total reads 60 million.

f1="${INPUT}"/SRR15852393_trimmed_1.fastq.gz
f2="${INPUT}"/SRR15852393_trimmed_2.fastq.gz
BASENAME=$(basename "${f1}" "_trimmed_1.fastq.gz")

# Run Salmon with GC bias and bootstrapping
salmon quant \
-i "${SALMON_INDEX}" \
-l A \
-1 "${f1}" \
-2 "${f2}" \
-p "${SALMON_THREADS}" \
--validateMappings --gcBias \
--numBootstraps 30 \
-o "${OUT}/${BASENAME}_salmon_quant"

elapsed=$SECONDS
echo "Completed salmon quantification for ${BASENAME}"
echo "Total elapsed time: $elapsed seconds" >> "${TIME_OUT}"
echo "Total elapsed time: $elapsed seconds"

```



#### Notes on alignment rates in Salmon and reasons 

The sample `SRR15852393` with **high GC-content and increased duplication levels and overepresented sequences** gad an alignment rate of 10% only. So salmon could not quasi-map mosty of the reads or I made some mistake in the trimming or indexing. I will re-check everything again.
- For paired-end sample of 2x28million reads with around 110bp length after trimming, and using decoy transcriptome index, gc-bias corrections, sequence bias correction, and bootsrapping, the alignment rate was only 12% and it took around 16 minutes with 7 threads. 
- It seems that the distribution of the read lengths affect the alignment rate as well. For that sample the distro was bimodal. If many 1/3 of reads are below 75bp after trimming, then using k=31 for the index could affect the alignment rate.
- Performing BLAST of the overepresented sequences give rRNA hits. Furtheremore 7 from 28 million reads were mapped to decoy sequences and discarded. This means that these reads were intronic, genomic, rRNA or other pre-mRNA sequences that could be mapped with higher confidence to the genome than the transcriptome.
- It is generally expected to have low alignment rates in Salmon compared to traditional aligners like STAR or HISAT2 that align to the genome and transciptome at the same time. Some report that Salmon alignment rates can be low. FPPE samples have reported to be of low quality. 
- The sample has veru high GC-content and overepresented sequences which could affect the alignment rate.
- The index was built with decoy transcriptome to account for spurious mappings. This reduced the alignment rate further.
- Maybe increasinng minimum read length and changing k-mer. The base quality is good so trimming parameters are not the issue. I use the cut_tail option in fastp to trim low-quality bases from the 3' end where quality usually drops. Maybe I should cut from both ends?
- For trimmed reads another solution is to reduce the `--minScoreFraction` from the default 0.65. It is discussed [here](https://github.com/COMBINE-lab/salmon/issues/533). No references yet for this. 


Output from Salmon
```
[2026-01-10 03:15:21.104] [jointLog] [info] Computed 189591 rich equivalence classes for further processing
[2026-01-10 03:15:21.104] [jointLog] [info] Counted 3380905 total reads in the equivalence classes
[2026-01-10 03:15:21.105] [jointLog] [info] Number of mappings discarded because of alignment score : 137506762
[2026-01-10 03:15:21.105] [jointLog] [info] Number of fragments entirely discarded because of alignment score : 7843579
[2026-01-10 03:15:21.105] [jointLog] [info] Number of fragments discarded because they are best-mapped to decoys : 7350523
[2026-01-10 03:15:21.105] [jointLog] [info] Number of fragments discarded because they have only dovetail (discordant) mappings to valid targets : 75102
```

I will test a sample `SRR15852395` with far worse quality and  `SRR15852396` with better length distribution and less duplication levels. The former was a nightmare, most reads did not meet the score threshold. The latter had an almost 30% alignemnt rate. Still many were mapped to decoys and discarded. 

```bash
# check the percent-mapping in all samples
grep -r "percent" ${HOME}/cancer_project/ngs_preprocess/5a_salmonout/

# check other info related to alignment 
grep -r -A 7 "num_processed" ${HOME}/cancer_project/ngs_preprocess/5a_salmonout/
```

some references to check for this issue:

### Mapping to genome with STAR aligner 

STAR always map to the genome and transcriptome simultaneously if the index was built with gtf annotation. It selects the best alignment for each read. STAR first aligns reads to entire genome, and only then searches for concordance between alignments and transcripts. This allows for discovery of splice variant in contrast to Salmon where every read is force to a coding region and much more strict constraints are applied to pseudo-alignemnt. <br>

Regardless of transcriptome, STAR projects the genome alignemnts to the transcriptome with the option `--quantMode TranscriptomeSAM`. The output is a BAM file with alignments projected to the transcriptome. This BAM file can be used for quantification with RSEM. One important question is what do you need the BAM file for? Depending on the downstream analysis, the BAM file must be sorted or unsorted. Here are some notes:

- For example, if I want to use Qualimap for QC of the alignments, the BAM file must be name sorted with samtools.
- For counting with featureCounts, the BAM file can be unsorted and should be genome-aligned (not transcriptome projected).
- For visualization in IGV, the BAM file must be sorted and indexed with samtools index.
- For variant calling with GATK, the BAM file must be sorted, indexed and marked for duplicates, then recalibration steps must follow.

Some Notes on the command:

- The `--runThreadN` option specifies the number of threads to use for parallel processing.
- The `--genomeDir` option specifies the path to the STAR genome index directory created
- The `--readFilesCommand zcat` option is used to read compressed FASTQ files directly.
- The `--outFileNamePrefix` option specifies the prefix for all output files generated by STAR for that sample. If the basename directory does not exist it will not be created automatically. For example, if the prefix is `path/to/output/basename_`, all the output files will start with that prefix. If you want a directory like `path/to/output/basename/basename_`, you need to create it `mkdir path/to/output/basename/` before running the STAR command.

```bash

f1=${HOME}/cancer_project/ngs_preprocess/2_trimmed/SRR15852393_trimmed_1.fastq.gz
f2="${f1/_1.fastq.gz/_2.fastq.gz}"
BASENAME=$(basename "${f1}" "_trimmed_1.fastq.gz")

# Example STAR mapping command for paired-end data
STAR --runThreadN 22 \
--genomeDir ${HOME}/cancer_project/ngs_preprocess/reference/hg38_star \
--readFilesIn ${f1} ${f2} \
--readFilesCommand zcat \ 
--outFileNamePrefix ${HOME}/cancer_project/ngs_preprocess/3_mapping/${BASENAME} \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts
```

The output files from STAR are:

- A per sample folder in the `3_mapping` directory with the name of the sample.
- Each sample folder contains multiple output files from STARand each folder contains the names of the sample as prefix. 
- Aligned.sortedByCoord.out.bam : BAM file sorted by coordinate with genome alignments.
- Aligned.toTranscriptome.out.bam : BAM file with alignments projected to the transcriptome for quantification.
- Log.final.out : summary of mapping statistics.

Informations about computer resources:

- For a sample with 30 million paired-end reads of 2x110bp length after trimming, STAR took around 4 minutes with 22 threads and used around 31GB of RAM. The index was built with splice junctions from GTF annotation.

### Post Alignment QC (RSeQC, Qualimap, MultiQC)

Starting with RSeQC. I will need the BAM file of each sample sorted by coordinate. Furthermore, I will need the gene model (.gtf file) in bed format. I have already done the conversion in the first section. Furthermore, I will need a subset of the genes for the gene body coverage function. The file I will uise is the one available from RSeQC.

```bash
# I should check the chromosome contigs names in BAM and BED files.
samtools idxstats ./3_mapping/SRR15852393/SRR15852393_Aligned.sortedByCoord.out.bam | head

head -n 5 ${HOME}/cancer_project/ngs_preprocess/reference/Homo_sapiens.GRCh38.115.bed

head -n 5 ${HOME}/cancer_project/ngs_preprocess/reference/Homo_sapiens.GRCh38.115.bed

```
An important note about using samtools to convert to `.bam.bai` files. When I run samtools and I used the .bam file path, the program for the RSeQC automatically found the `.bai` when it was in the same directory. 

The time required to run the RSeQC for one sample is around 12 minutes. Qualimap requires the bam file to sorted buy coordinate and name (it automatically sorts by read name if i choose the -pe option). However, if the bam is not sorted by name, Qualimap's single-threaded mode is very slow, taking around 1200 minutes. <br>
Instead it would be better to create a **temporary bam_name_sort** file before we run qualimap. A bam file PE with 30 million reads took 116 seconds to be sorted by name and the bam file has a size of 1.3GB. The .bai version is around 7MB and it's coordinate-sorted version around 6GB. The transcriptome-Sam file around 1.5GB. The qualimap took 477 seconds with the name-sorted bam file as input. <br>

Below I have kept a part of my posty_alignqc bash script that runs Qualimap for each sample to remember the details. 

```bash
#!/usr/bin/env bash 
MAP_DIR="${HOME}"/cancer_project/ngs_preprocess/3_mapping
SALMON_DIR="${HOME}"/cancer_project/ngs_preprocess/5a_salmonout
QC_DIR="${HOME}"/cancer_project/ngs_preprocess/4_qualimap
REF_DIR="${HOME}"/cancer_project/ngs_preprocess/reference

# necessary files for RSeQC
BED_HOUSEGENES="${REF_DIR}/hg38.HouseKeepingGenes.bed"
BED_ANNOTATION="${REF_DIR}/Homo_sapiens.GRCh38.115.bed"

samtools sort -n -@ 16 -o "${QC_DIR}/${sample_name}.namesort.bam" "$bam"

qualimap rnaseq \
  -bam "$bam" \
  -gtf "${REF_DIR}/Homo_sapiens.GRCh38.115.gtf" \
  -outdir "${QC_DIR}/${sample_name}" \
  -outformat HTML \
  --java-mem-size=16G >& "${QC_DIR}/${sample_name}/${sample_name}_qualimap.log"

# remove temporary name_sorted bam
rm "${QC_DIR}/${sample_name}.namesort.bam"

```

**Verdict**
Since RSeQC takes time but not too much RAM or threads, it would be possible to parallelize the analysis of specific samples, around 4-5 with possibly 5 cores each. Qualimap will run with temporary name-sort bam files either parallelized by 2 or in a for loop (depending how RAM handles both Samples at the same time). 

### Quantification with RSEM 
RSEM requires a reference index built with the same genome fasta and gtf annotation file as STAR (laready done). RSEM can take as input either the FASTQ files directly or the BAM file from STAR with transcriptome alignments. Here I will use the latter approach. 
- An important note [here](https://deweylab.biostat.wisc.edu/rsem/rsem-calculate-expression.html?utm_source=chatgpt.com) is that RSEM requires the BAM file to be sorted by read name when using STAR transcriptome alignments as input. This is because RSEM expects paired-end reads to be adjacent in the BAM file for proper quantification. If the BAM file is sorted by genomic coordinate, RSEM may not correctly identify read pairs, leading to inaccurate expression estimates. (MIGHT BE INCORRECT)
- With 12 cores for ~30 million sam-tools took 31 seconds to sort the BAM file by read name.
- I noticed that when I used the `rsem-sam-validator` tool, the original transcriptome BAM file from STAR was valid for RSEM and the sorted gave an error `map1 read name not equal to map2 read name`. I run the unsorted Transcriptome BAM file with RSEM and it took 272 seconds with 11 threads for quantification. I should perform correlation analysis between Salmon and RSEM at the gene level and make a scatterplot with a linear fit and R^2 value.

```bash
# check with samtools, output: SO:unsorted or not at all
samtools view -H Aligned.toTranscriptome.out.bam | grep '^@HD'

# the command to sort by read name 
MAP_DIR="${HOME}/cancer_project/ngs_preprocess/3_mapping"
RSEM_REF="${HOME}/cancer_project/ngs_preprocess/reference/hg38_rsem"
OUT_DIR="${HOME}/cancer_project/ngs_preprocess/5b_starout"

SECONDS=0
samtools sort -n -@ 12 -o "${MAP_DIR}/SRR15852393/SRR15852393_transcriptome_name_sorted.bam" \
"${MAP_DIR}/SRR15852393/SRR15852393_Aligned.toTranscriptome.out.bam"
elapsed=$SECONDS
echo "Sorted BAM by read name in $elapsed seconds"

# sanity checks
samtools view -H "${MAP_DIR}/SRR15852393/SRR15852393_transcriptome_name_sorted.bam" | grep '^@HD'

# the input was invalud here 
rsem-sam-validator "${MAP_DIR}/SRR15852393/SRR15852393_transcriptome_name_sorted.bam"

# the input is valid here
rsem-sam-validator "${MAP_DIR}/SRR15852393/SRR15852393_Aligned.toTranscriptome.out.bam" 

# if it does not pass the validation, lets try this:
samtools collate -@ 12 -o "${MAP_DIR}/SRR15852393/SRR15852393_collated.bam" \
"${MAP_DIR}/SRR15852393/SRR15852393_transcriptome_name_sorted.bam"

# RSEM quantification command and time it 
SECONDS=0
mkdir -p "${OUT_DIR}/SRR15852393_rsem"
rsem-calculate-expression --paired-end --alignments \
--strandedness reverse \
-p 11 \
"${MAP_DIR}/SRR15852393/SRR15852393_Aligned.toTranscriptome.out.bam" \
"${RSEM_REF}/hg38_rsem" \
"${OUT_DIR}/SRR15852393_rsem/SRR15852393_rsem"
elapsed=$SECONDS
echo "Completed RSEM quantification for SRR15852393 in $elapsed seconds"

```
--outSAMattributes NH HI AS nM


## Some technical details in GNU parallel
An important note for GNU parallel is to understand where is the job template and what is part of a specific option. The job template could include several commands separated by semicolon `;`. The job template is inside single quotes `'...'`. The input values per job can be represented by placeholder variables such as `{}` or `{1}`, `{2}` etc. depending on the number of input sources. The input values are provided after the `:::` operator or by a pipe if there is one preceding. <br>

Other useful options are:
- The `--jobs` or `-j` option specifies the maximum number of jobs to run in parallel. For example, `-j 4` will run up to 4 jobs simultaneously.
- The `--tag` option adds the name of the input file before prepending the output of the command. An example was in prefetch SRAs, where in the latter case you will get the usual output plus the SRR ID.
- The `--linebuffer` option ensures that the output from each job is printed in its entirety before the output from the next job begins. This prevents interleaving of output from different jobs, which can make it difficult to read.
- The `--verbose` option provides detailed information about the execution of the jobs, including when each job starts and finishes. This can be useful for debugging and monitoring the progress of the parallel tasks.
- The `--tagstring` option allows you to customize the tag that is prepended to the output of each job. You can include placeholders such as `{#}` for the job number and `{1}`, `{2}`, etc., for the input values. This is useful for identifying which output corresponds to which input.

Some examples of how parallel works:

```bash

parallel echo ::: A B C D ::: 1 2 3

# this will pair the inputs 
parallel --link echo ::: A B C D ::: 1 2 3 4

# using placeholders variables
parallel echo {1} {2} ::: A B C D ::: 1 2 3
parallel --linebuffer --tagstring "job {#} with ({1} {2})" echo ::: A B C D ::: 1 2 3 

# use multiple commands per job, comma for command sep and bash -c
parallel --linebuffer --tagstring "job {#} with ({1} {2})" \
    bash -c 'echo "Processing {1} and {2}"; echo "Ending with {1}_{2}"' \
    ::: A B C D ::: 1 2 3 4 5

# put a maximum number of jobs to run in parallel the above
# the bit in the single quotes is the job template 
parallel --jobs 10 --linebuffer \
  --tagstring 'job {#} with ({1} {2})' \
  'echo "Processing {1} and {2}"; echo "Ending with {1}_{2}"' \
  ::: A B C D ::: 1 2 3
```

A more specific example with a custom script which has a duration of 2 seconds and placeholder variables. In --verbose, we can observe the moment that a job finishes another starts.

```bash custom.sh
#!/usr/bin/env bash

JOB_ID=${1} # this is placeholder
sleep 2
echo "Finished the job #${JOB_ID}"
```
```bash custom2.sh
#!/usr/bin/env bash

sleep ${1} # placeholder for first argument 
echo "The object ${2} and slept for ${1} seconds"

```

```bash 
# works 5 jobs each time, placeholder gets the value from input 
parallel -j 5 ./custom.sh {} ::: 1 2 3 4 5 6 7

# here no placeholder is used, so the input values are passed as arguments
parallel --jobs 5 --verbose ./custom.sh ::: 1 2 3 4 5 6 7 8 9 10

parallel --jobs 5 --verbose ./custom.sh ::: $(seq 10) # same as above

# get job from custom2 and tagstring too
parallel -j 4 --linebuffer --tagstring "job {#}" \
./custom2.sh ::: $(seq 1 3) ::: A B C D

# specify arguments irrespective of the number of input sources, if j =1 for loop
parallel -j 12 --linebuffer --tagstring "job {#} from ({1} {2})" \
./custom2.sh {2} {1} ::: A B C D ::: 1 2 3

```

When I am getting an STDIN from pipe operator I am using the `{}` placeholder to represent each line of input. No numbers are used like the bash scripts where I specify the inputs. For example:

```bash
ls sra_data/parallel_fetch/ | cat | parallel -j 4 --tagstring "{#} working with {}" \ echo "actually using the {}"
```
## Advanced uses of Commands 

### Tricks with sed and regex
- In this example I am using find sort and sed to extract pattern or remove pattern from filenames.
- Imagine if you want to remove some pattern from the fastq file and want to extract a basename. The example below is for `./1_fastq/SRR15852393`.


```bash
find ./1_fastq/ -iname "*.zip" | sed -E "s/_[12]_fastqc\.zip$//" | sed -E "s#\./[1-9]_fastq/##"

# I can also combine the two sed commands with ; 
find ./1_fastq/ -iname "*.zip" | sed -E "s/_[12]_fastqc\.zip$//; s#\./[1-9]_fastq/##"

```
What happned here?
- The first sed command substitutes `s/REGEX/SUBSTITUTE/END` with nothing as substitute. In the second sommand how everything works:
    - `s#\./[1-9]_fastq/##` : the `#` is used instead of `/` to avoid confusion with the slashes in the regex. The regex `\./[1-9]_fastq/` matches the pattern starting with `./`, followed by a digit from 1 to 9, then `_fastq/`. The backslash before the dot `\.` is used to escape the dot, indicating that we are looking for a literal dot character rather than any character.
    - Why do i not escape the `/` ? Because in this case `/` is not a special character in regex, so it does not need to be escaped. The regex engine interprets it as a literal forward slash.
    - Why am i using two `##` ? The first `#` is the delimiter for the `s` command, and the second `#` indicates the end of the regex pattern. The final `#` indicates the end of the substitute part, which is empty in this case. It does not mean thhe end of the line whose patter is being matched. 


### Awk for text processing

```bash 
awk -F '[[:space:]]*\\|[[:space:]]*' 'BEGIN { OFS = "\t"; print "method", "cores", "time" } {print $2, $3, $4}' ./scripts/exploratory/time.log | sed 's#seconds##g' > results.txt

```

difference in basename uses:
In the second case the dirname returns the parent directory of the file path, which is `SRR123` in this case. Then basename extracts the last part of that path, which is again `SRR123`. So both approaches give the same result here.

```bash

# hypothetical directory with files
bam_dir="/data/3_mapping/SRR123/SRR123_Aligned.sortedByCoord.out.bam"

base=$(basename "${bam_dir}" "_Aligned.sortedByCoord.out.bam")
echo "${base}"  # Output: SRR123

base2=$(basename $(dirname "${bam_dir}"))
echo "${base2}"  # Output: SRR123)

```
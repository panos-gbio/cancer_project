#!/usr/bin/env bash

# ADD PARAMETERS AND APPROPRIATE PATHS 

JOBS=10
SRA_FILE=~/practice/sra_data/SRR_Acc_List.txt
OUTPUT_FETCH=~/practice/sra_data/parallel_fetch
OUTPUT_FASTQ=~/practice/sra_data/parallel_fastq

# functions from SRA toolkit 
PREFETCH=~/Downloads/sratoolkit.3.2.1-ubuntu64/bin/prefetch
FASTERQ_DUMP=~/Downloads/sratoolkit.3.2.1-ubuntu64/bin/fasterq-dump


# Check that SRA list file exists and is not empty
if [[ ! -s "${SRA_FILE}" ]]; then
    echo "SRA text file is empty or not found."
    exit 1
else
    echo "SRA text file list exists and there are $(cat "${SRA_FILE}" | wc -l) SRAs to fetch."
fi 


# check that the output directories exist
if [[  ! -d "${OUTPUT_FETCH}" || ! -d "${OUTPUT_FASTQ}" ]]; then
    echo "You have not specified the correct output directory for prefetch or fastq-dump."
    exit 1
elif [[ -d "${OUTPUT_FETCH}" && -d "${OUTPUT_FASTQ}" ]]; then
    echo "Output directories for prefetch and fastq-dump exist."
fi

# check out the content of the SRA list file
while IFS= read -r read; do
    echo "Fetching ${read}..."
done < "${SRA_FILE}" # redirecting variable expansion 


# Parallelize prefetching with GNU parallel

# cat retains line breaks and GNU recognises each line as a separate input job
if [[ -z "$(ls -A "${OUTPUT_FETCH}")" ]]; then # check is the output dir is empty 
    echo "Output directory is empty, starting fetch..."
    cat "${SRA_FILE}" \
        | parallel -j "${JOBS}" --linebuffer --tag \
            "${PREFETCH}" -O "${OUTPUT_FETCH}" --max-size u {}
    echo "Fetching completed."

else
    echo "Output directory is not empty, no need to fetch. Remove anything but SRR directories."
fi

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
    echo "FASTQ output directory is not empty, no need to run fasterq-dump. Remove all anything non FASTQ."
fi

# either use GNU parallel or a tool to speed up the zipping of fastq files
#TODO: add gzip parallelization
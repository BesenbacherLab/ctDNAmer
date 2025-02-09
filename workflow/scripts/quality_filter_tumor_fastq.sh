#!/usr/bin/env bash

exec 2> "${snakemake_log[0]}" # redirect errors
echo "errors redirected"  >> "${snakemake_log[0]}" 

input_files=(${snakemake_input[input_files]}) # read in text file with input file names
echo "input file with the file list: $input_files" >> "${snakemake_log[0]}" 

while read -r line || [[ -n "$line" ]]; do # loop over input files
    filename="$line"
    echo "working on file: $filename"  >> "${snakemake_log[0]}"
    outfile=${filename%.*}
    echo "writing to filtered file: $outfile.filtered.fastq"  >> "${snakemake_log[0]}"
    seqtk trimfq -b 5 -e 5 $filename | seqtk seq -q37 -n N > "$outfile.filtered.fastq" # remove 5 bp from each end and set all bases with BQ < 37 to N
    echo "$outfile.filtered.fastq" >> "${snakemake_output[0]}"
done < "$input_files"
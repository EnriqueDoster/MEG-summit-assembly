#!/usr/bin/env bash
# Run this script to make the various SBATCH files needed for al your samples
# Use: parallel -j 1 "../parallel_assembly.bash {}" ::: input_dir/*R1.fastq.gz
out_directory="assembly"
input_dir

samplename=$( basename $1 | sed -r 's/R1_001.fastq.gz//' )
input_dir=$(dirname $1)

echo "${samplename}"

python bin/sbatch_assembly_master.py -f ${input_dir}${samplename}R1_001.fastq.gz -r ${input_dir}${samplename}R2_001.fastq.gz -o ${out_directory}

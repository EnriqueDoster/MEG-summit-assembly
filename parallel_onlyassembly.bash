#!/usr/bin/env bash
# Run this script to make the various SBATCH files needed for each sample and output to another directory
# Usage: parallel -j 1 "PATH/TO/MEG-summit-assembly/parallel_onlyassembly.bash {}" ::: /PATH/TO/READS/*R1.fastq.gz

out_directory="assembly" # make this directory before running this script, no "/" at the end

samplename=$( basename $1 | sed -r 's/R1.fastq.gz//' )
input_dir=$(dirname $1)

echo "${samplename}"

python bin/sbatch_onlyassembly_master.py -f ${input_dir}${samplename}R1.fastq.gz -r ${input_dir}${samplename}R2.fastq.gz -o ${out_directory}

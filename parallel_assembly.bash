#!/usr/bin/env bash
# Run the script from the location of the sam files and script must be in a different directory
# Use: parallel -j 1 "../parallel_assembly.bash {}" ::: ~/Dropbox/Reference_resources/Useful_code/Python/*R1.fastq.gz
out_directory="assembly"

samplename=$( basename $1 | sed -r 's/R1.fastq.gz//' )

echo "${samplename}"

python sbatch_assembly_master.py -f ${samplename}R1.fastq.gz -r ${samplename}R2.fastq.gz -o ${out_directory}

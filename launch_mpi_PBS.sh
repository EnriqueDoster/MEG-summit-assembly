#!/bin/bash
#qsub -N amrplusplus
#qsub -q small
#qsub -l procs=1
#qsub -l walltime 20:00:00
#qsub -l M enriquedoster@gmail.com

module purge
module load java/jdk1.8.0_144
module load singularity/current
module load ompi/gnu.mesabi

mpirun --pernode ./nextflow run main_amrplusplus_snp.nf -resume -profile slurm \
-w /home/noyes046/slizo001/projects/cipro/work --threads 15 \
--output /home/noyes046/slizo001/projects/cipro/cipro_amrplusplus --host /home/noyes046/slizo001/projects/HostGenome/GRCh38_latest_genomic.fna.gz \
--reads "/home/noyes046/slizo001/projects/cipro/data/*_{1,2}.fastq.gz" -with-mpi

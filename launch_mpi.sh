#!/bin/bash
#SBATCH --job-name=nf_HMM
#SBATCH --partition=shas
#SBATCH --ntasks=1
#SBATCH --qos=long
#SBATCH --cpus-per-task=1
#SBATCH --time=100:00:00
#SBATCH --export=ALL
#SBATCH --mail-user=enriquedoster@gmail.com
#SBATCH --mail-type=ALL

module purge
module load jdk/1.8.0
module load singularity/2.5.2
module load gnu_parallel/20160622
module load gcc/6.1.0
module load openmpi/2.0.1

mpirun --pernode ./nextflow run main_HMM_SNP_dedup.nf -resume -profile slurm \
-w /scratch/summit/edoster@colostate.edu/FCmeat_results/work --threads 15 \
--output /scratch/summit/edoster@colostate.edu/FCmeat_results --host /scratch/summit/edoster@colostate.edu/Livestock_complete_genomes.fa \
--reads "/scratch/summit/edoster@colostate.edu/FC_meat/*_R{1,2}.fastq.gz" -with-mpi

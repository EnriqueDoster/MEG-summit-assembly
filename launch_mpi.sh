#!/bin/bash
#SBATCH --job-name=nf_HMM
#SBATCH --partition=smem
#SBATCH --ntasks=1
#SBATCH --qos=long
#SBATCH --cpus-per-task=1
#SBATCH --time=50:00:00
#SBATCH --export=ALL
#SBATCH --mail-user=enriquedoster@gmail.com
#SBATCH --mail-type=ALL

mpirun --pernode ./nextflow run main_HMM_SNP_dedup.nf -resume -profile slurm -w /scratch/summit/edoster@colostate.edu/SNP_HMM/proj7/work --threads 1 --output /scratch/summit/edoster@colostate.edu/SNP_HMM/proj7/ --host /scratch/summit/edoster@colostate.edu/HMM_steps/MEG-summit-assembly/Livestock_complete_genomes.fa --reads "/scratch/summit/edoster@colostate.edu/proj7/AFE1*_R{1,2}.fastq.gz" -with-mpi 


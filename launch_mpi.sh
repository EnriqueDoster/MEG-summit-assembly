#!/bin/bash
#SBATCH --partition=shas
#SBATCH --ntasks=1
#SBATCH --qos=normal
#SBATCH --cpus-per-task=1
#SBATCH --time=20:00:00
#SBATCH --export=ALL
#SBATCH --mail-user=enriquedoster@gmail.com
#SBATCH --mail-type=ALL

mpirun --pernode ./nextflow run main_dedup.nf -w /scratch/summit/edoster@colostate.edu/HMM_steps/proj7/work --threads 1 --output /scratch/summit/edoster@colostate.edu/HMM_steps/proj7/ --host /scratch/summit/edoster@colostate.edu/HMM_steps/MEG-summit-assembly/Livestock_complete_genomes.fa --reads "/scratch/summit/edoster@colostate.edu/proj7/AFE13_R{1,2}.fastq.gz" -profile slurm -with-mpi -resume
#srun ./nextflow run main_dedup.nf -w /scratch/summit/edoster@colostate.edu/HMM_steps/proj7/worksrun --threads 1 --output /scratch/summit/edoster@colostate.edu/HMM_steps/proj7_srun/ --host /scratch/summit/edoster@colostate.edu/HMM_steps/MEG-summit-assembly/Livestock_complete_genomes.fa --reads "/scratch/summit/edoster@colostate.edu/proj7/AFE13_R{1,2}.fastq.gz" -profile slurm -with-mpi -resume


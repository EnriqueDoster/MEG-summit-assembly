# MEG-summit-assembly
### Written by:	 Enrique Doster and Steven Lakin
### Date:		     2 October 2018
### Purpose: 	   These set of repository is for making sbatch scripts for each pair of fastq samples and submits those scripts to the Slurm job scheduler

### The assembly process can begin with either raw reads or non-host reads.
- Carefully consider which scripts are being used and confirm that the right number of nucleotides are removed from the start of the reads. 
- First, ensure that you samples reads are uploaded to SUMMIT, in /scratch/summit/$USER/ . I recommend logging in to ANGUS and using sftp to "put" reads unto SUMMIT. CSU's VPN makes it difficult to sftp from SUMMIT to ANGUS.

# MEG-summit-assembly Workflow
### You will have to edit the parallel_onlyassembly.bash and SLURM_task_launcher.sh scripts
### Login to SUMMIT
ssh user@colostate.edu@login.rc.colorado.edu

## To access Summit compile nodes:
ssh scompile
## Purge and load existing modules
module purge

module load jdk/1.8.0

module load singularity/2.5.2

module load gnu_parallel/20160622

module load gcc/6.1.0

module load openmpi/2.0.1

## Navigate to directory where you want stage your scrips
cd /scratch/summit/$USER/
## install nextflow
curl -s https://get.nextflow.io | bash
## Download this repository
git clone https://github.com/EnriqueDoster/MEG-summit-assembly.git
cd MEG-summit-assembly


### View queue
squeue -u $USER
### Cancel all the jobs for a user:
scancel -u $USER


### Still in progress - development of nextflow script
# main_dedup.nf
# use the following command in a SLURM launcher script.
mpirun --pernode ./nextflow run main_dedup.nf -w /scratch/summit/edoster@colostate.edu/HMM_steps/proj7/work --threads 1 --output /scratch/summit/edoster@colostate.edu/HMM_steps/proj7/ --host /scratch/summit/edoster@colostate.edu/HMM_steps/MEG-summit-assembly/Livestock_complete_genomes.fa --reads "/scratch/summit/edoster@colostate.edu/proj7/AFE13_R{1,2}.fastq.gz" -profile slurm -with-mpi -resume

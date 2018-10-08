# MEG-summig-assembly
### Written by:	 Enrique Doster and Steven Lakin
### Date:		     2 October 2018
### Purpose: 	   These set of repository is for making sbatch scripts for each pair of fastq samples and submits those scripts to the Slurm job scheduler

### The assembly process can begin with either raw reads or non-host reads. Carefully consider which scripts are being used and confirm that the right number of nucleotides are removed from the start of the reads. 

### First, ensure that you samples reads are uploaded to SUMMIT, in /scratch/summit/$USER/ . I recommend logging in to ANGUS and using sftp to "put" reads unto SUMMIT. CSU's VPN makes it difficult to sftp from SUMMIT to ANGUS.

# MEG-summit-assembly Workflow
### You will have to edit the parallel_onlyassembly.bash and SLURM_task_launcher.sh scripts
### Login to SUMMIT
ssh user@colostate.edu@login.rc.colorado.edu

## To access Summit compile nodes:
ssh scompile

## Purge and load existing modules
module purge
module load gnu_parallel/20160622

## Navigate to directory where you want stage your scrips
cd /scratch/summit/$USER/
git clone https://github.com/EnriqueDoster/MEG-summit-assembly.git
cd MEG-summit-assembly

## Edit MEG-summit-assembly/parallel_onlyassembly.bash file to change the output_directory and confirm that your samples have similar naming conventions (ie. sample_R1_001.fastq.gz vs sample_R1.fastq.gz)
### Create sbatch scripts for your samples
parallel -j 1 "/scratch/summit/$USER/MEG-summit-assembly/parallel_onlyassembly.bash {}" ::: /scratch/summit/$USER/proj_dir/*R1.fastq.gz

## Test one sbatch script in your output_directory
#View queue
squeue -u $USER
#Cancel all the jobs for a user:
scancel -u $USER
#If sbatch script completes succesfully, prepare the SLURM_task_launcher.sh script. Edit output dir and timing of submission

## Submit the remaining sbatch scripts using the SLURM_task_launcher.sh
sbatch /scratch/summit/$USER/MEG-summit-assembly/SLURM_task_launcher.sh

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

## Navigate to directory where you want stage your scrips
cd /scratch/summit/$USER/

## Download this repository, unzip HMM models, and change permissions in directory
git clone https://github.com/EnriqueDoster/MEG-summit-assembly.git

cd MEG-summit-assembly

gunzip containers/data/HMM/*

chmod 755 -R *

## Outline of pipeline installation and execution:

- Java is necessary to install nextflow and run the pipeline
usually this can be done with by loading modules available to all users on a server

  - $ module load jdk/1.8.0

- Other modules required, but will be loaded automatically:

  - Singularity - contains all bioinformatic software required for running the pipeline- https://github.com/singularityhub/singularityhub.github.io/wiki \http://singularity-hub.org/ | 
  
    - $ module load singularity/2.5.2
 
  - GCC compiler

    - $ module load gcc/6.1.0

  - OpenMPI

    - $ module load openmpi/2.0.1

  - Once the repository is downloaded, I recommend going into that directory and installing nextflow (with java loaded)

    - $ git clone https://github.com/EnriqueDoster/MEG-summit-assembly.git

    - $ cd MEG-summit-assembly/

    - $ curl -s https://get.nextflow.io | bash

  - Next, two housekeeping commands to unzip files and make executable

    - $ gunzip containers/data/HMM/*

    - $ chmod 755 -R *

 - The "launch_mpi.sh" script must be edited.
 
   - Edit #SBATCH options

   - Edit last command on file to specify the following flags:
 
     - --reads | path to input fastq files

     - --host | path to host genome fasta

     - --output | output path

     - --threads | number of threads to use (cores in a node)

     - -w | path for work directory

     - -profile | server configuration profiles

 - When choosing the locations for the output, please consider that the output directory and work directory take up the same amount of space and can overwhelm disk quotas. 

   - For a metagenomic sample with paired reads (7GB), we estimate 93GB of output to be stored and 93GB of temporary files (total 186 GB/sample).

   - Once the whole pipeline is complete, the work directory can be deleted. 

   - We recommend running the pipeline from within the MEG-summit-assembly/ directory, but the --output and -w directories can be in distinct locations on the server. 

   - /scratch/ directories on shared servers are typically purged after a certain time period (~90 days) but are typically the best option for the work directory. Check with IT support for further details.

 - Edit the "nextflow.config" file
 
   - Review/edit quality control parameters
 
 - Review the configuration file and the "maxForks" variable to control how many samples will be run in parallel

   - $ nano  config/slurm.config

   - It's important to change the name of the partitions and queue names for your server.

 - The "launch_mpi.sh" script is submitted using sbatch.

   - $ sbatch launch_mpi.sh

 - The launch_mpi.sh script is submitted as a single process (slurm job) and once picked up from the queue will organize the submission of further jobs to complete the pipeline. 

   - maxForks controls how many processes (ie slurm jobs) will be submitted at a time.

   - This must be balanced with the use of the "--threads" flag to optimize the run time. 
 
   - Depending on the duration of the run, the main process for the launch_mpi.sh may be dropped by the server

    - Fortunately, nextflow allows the use of the "-resume" flag (included by default) to pick up where it left off. 

    - $ sbatch launch_mpi.sh
 
 - Run times are highly variable, but on a server running 24 threads it could take up 14 hours per sample (3GB zipped paired fastq) when running the longest pipeline "main_HMM_SNP_dedup.nf". The "maxForks" flag in the config/slurm.config file is your friend. The "main_HMM_base.nf" is the most streamlined version that only includes the amrplusplus version with HMMs (instead of bwa alignment).

   - Once the run is completed, you should get an email and can verify by viewing the ".nextflow.log" file

   - You can now delete the work "-w"  directory to save space. 

- Basic troubleshooting

   - If running on slurm, some individual steps of the pipeline will be unable to complete due to memory or time constraints. You can edit the "config/slurm.config" file to change this.
    - For example, the AssembleReads step in the pipeline can take more memory than expected. The SBATCH flag "--ntasks-per-node" can be modified to request more cores of RAM for that particular process. Some processes reu

    - Then, you can simply submit the launch_mpi.sh script again.


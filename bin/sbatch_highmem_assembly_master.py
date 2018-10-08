## Python script that receives a sample's forward and reverse read locations, creates SBATCH script

## Import packages
import argparse
import os
import sys
import gzip
import io

#USAGE : python sbatch_assembly_master.py -f ~/Dropbox/Reference_resources/Useful_code/Python/test_R1.gz -r ~/Dropbox/Reference_resources/Useful_code/Python/test_R2.gz -o test/

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--forward', type=str, required=True,
    	                help='Forward input file')
parser.add_argument('-r', '--reverse', type=str, required=True,
    	                help='Reverse input file')
parser.add_argument('-o', '--output', required=True,
    	                help='Specify output directory for SBATCH script files')


## Other parameters
TRIMMOMATIC = '/opt/Trimmomatic-0.36'
adapters = "../containers/data/adapters/nextera.fa"
fqc_adapters = "../containers/data/adapters/nextera.tab"
leading = '10' ## To remove UMI's
trailing = '3'
slidingwindow = "4:15"
minlen = '36'
threads= '8'



if __name__ == "__main__":
    args = parser.parse_args()
    if not os.path.exists(args.output):
        print("outdir path doesn't exist. trying to make")
        os.makedirs(args.output)
    print(args.forward.split('/')[-1].split('_')[0])
    samplename = str(args.forward.split('/')[-1].split('_R1')[0])
    print(samplename)
    with open(args.output + '/' + samplename + '_SLURM_script.sh', 'w') as fout:
        ## First, write the SBATCH information
        shebang = '#!/bin/bash'
        job_name = ('#SBATCH --job-name={}_assembly%j'.format(samplename)) 
        partition = '#SBATCH --partition=smem'
        nodes= '#SBATCH --nodes=1'
	mem = '#SBATCH --mem=60000'
	log = ('#SBATCH --output={}_log%j'.format(samplename))
        time = '#SBATCH --time=23:55:00'
        fout.write('{}\n{}\n{}\n{}\n{}\n{}\n{}\n'.format(shebang,job_name,partition,nodes,mem,log,time)) # Write out all necessary information for slurm at top of script
        # Load modules
        purge_mod = 'module purge'
        mod_jkd = 'module load jdk/1.8.0'
        mod_singularity = 'module load singularity/2.5.2'
        mod_gcc = 'module load gcc/6.1.0'
        mod_mpi = 'module load openmpi/2.0.1'
        #execute mod_singularity
        #exec_singularity = 'srun singularity run /scratch/summit/edoster@colostate.edu/EnriqueDoster-MEG-summit-assembly-master-latest.simg'
        # Change to output directory
        cd_out_dir= ('cd {}'.format(args.output))
        fout.write('{}\n{}\n{}\n{}\n{}\n{}\n'.format(purge_mod,mod_jkd,mod_singularity,mod_gcc, mod_mpi, cd_out_dir)) # Write out all necessary information for slurm at top of script
            ## QC control
        trimm = ('singularity exec /scratch/summit/edoster@colostate.edu/EnriqueDoster-MEG-summit-assembly-master-latest.simg /usr/lib/jvm/java-7-openjdk-amd64/bin/java -jar {}/trimmomatic-0.36.jar \
                PE \
                {} {} -baseout {} \
                -threads {} \
                ILLUMINACLIP:{}:2:30:10:3:TRUE \
                LEADING:{} \
                TRAILING:{} \
                SLIDINGWINDOW:{} \
                MINLEN:{} \
                2> {}.trimmomatic.stats.log'.format(TRIMMOMATIC,args.forward,args.reverse,samplename,threads,adapters,leading,trailing,slidingwindow,minlen,samplename))
        clean_up_trimm = ('gzip -c {}_1P > {}.1P.fastq.gz \n \
            gzip -c {}_2P > {}.2P.fastq.gz \n \
            rm {}_1P \n \
            rm {}_2P \n \
            rm {}_1U \n \
            rm {}_2U'.format(samplename,samplename,samplename,samplename,samplename,samplename,samplename,samplename))
        fout.write('{}\n{}\n'.format(trimm,clean_up_trimm))
        ## Assembly commands
        mk_temp_idba = ('mkdir -p temp{}/idba'.format(samplename))
        merge_fq = ('singularity exec /scratch/summit/edoster@colostate.edu/EnriqueDoster-MEG-summit-assembly-master-latest.simg /usr/local/bin/fq2fa --merge --filter <( zcat {}.1P.fastq.gz ) <( zcat {}.2P.fastq.gz ) temp{}/interleavened.fasta'.format(samplename,samplename,samplename))
        assemble_fa = ('singularity exec /scratch/summit/edoster@colostate.edu/EnriqueDoster-MEG-summit-assembly-master-latest.simg /usr/local/bin/idba_ud --num_threads {} -r temp{}/interleavened.fasta -o temp{}/idba'.format(threads,samplename,samplename))
        rename_contig = ('cp temp{}/idba/contig.fa {}.contigs.fasta'.format(samplename,samplename))
        fout.write('{}\n{}\n{}\n{}\n'.format(mk_temp_idba,merge_fq,assemble_fa,rename_contig))


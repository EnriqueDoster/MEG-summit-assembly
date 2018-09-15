## Python script that receives a sample's forward and reverse read locations, creates SBATCH script

## Import packages
import argparse
import os
import sys
import gzip
import io


## Input = 2 samples

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
adapters = "containers/data/adapters/nextera.fa"
fqc_adapters = "containers/data/adapters/nextera.tab"
leading = '3'
trailing = '3'
slidingwindow = "4:15"
minlen = '36'
threads= '22'



if __name__ == "__main__":
    args = parser.parse_args()
    if not os.path.exists(args.output):
        print("outdir path doesn't exist. trying to make")
        os.makedirs(args.output)
    print(args.forward.split('/')[-1].split('_')[0])
    samplename = str(args.forward.split('/')[-1].split('_')[0])
    print(samplename)
    with open(args.output + '/' + samplename + '_SLURM_script.sh', 'w') as fout:
        ## First, write the SBATCH information
        shebang = '#!/bin/bash'
        line2 = '#SBATCH --job-name=launcher'
        line3 = '#SBATCH --output=assembly_log%j'
        line4 = '#SBATCH --cpus=1'
        line5 = '#SBATCH --ntasks-per-node=20'
        line6 = '#SBATCH --time=20:00:00'
        fout.write('{}\n{}_{}\n{}\n{}\n{}\n{}\n'.format(shebang,line2,samplename,line3,line4,line5,line6)) # Write out all necessary information for slurm at top of script
        # Load modules
        purge_mod = 'module purge'
        mod_jkd = 'module load jdk/1.8.0'
        mod_singularity = 'module load singularity/2.5.2'
        #execute mod_singularity
        exec_singularity = 'singularity run /projects/edoster@colostate.edu/meg-assembly/summit-assembly/lakinsm-summit-assembly-latest.img'
        # Change to output directory
        cd_out_dir= ('cd {}'.format(args.output))
        fout.write('{}\n{}\n{}\n{}\n{}\n'.format(purge_mod,mod_jkd,mod_singularity,exec_singularity,cd_out_dir)) # Write out all necessary information for slurm at top of script
        ## QC control
        trimm = ('/usr/lib/jvm/java-7-openjdk-amd64/bin/java -jar {}/trimmomatic-0.36.jar \
                PE \
                -threads {} \
                {} {} -baseout {} \
                ILLUMINACLIP:{}:2:30:10:3:TRUE \
                LEADING:{} \
                TRAILING:{} \
                SLIDINGWINDOW:{} \
                MINLEN:{} \
                2> {}.trimmomatic.stats.log'.format(TRIMMOMATIC,threads,args.forward,args.reverse,samplename,adapters,leading,trailing,slidingwindow,minlen,samplename))
        clean_up_trimm = ('gzip -c {}_1P > {}.1P.fastq.gz \n \
            gzip -c {}_2P > {}.2P.fastq.gz \n \
            rm {}_1P \n \
            rm {}_2P \n \
            rm {}_1U \n \
            rm {}_2U'.format(samplename,samplename,samplename,samplename,samplename,samplename,samplename,samplename))
        fout.write('{}\n{}\n'.format(trimm, clean_up_trimm))
        ## Assembly commands
        mk_temp_idba = ('mkdir -p temp{}/idba'.format(samplename))
        merge_fq = ('fq2fa --merge --filter  <( zcat {}.1P.fastq.gz) <( zcat {}.2P.fastq.gz ) temp{}/interleavened.fasta'.format(samplename,samplename,samplename))
        assemble_fa = ('idba_ud --num_threads {} -l temp{}/interleavened.fasta -o temp{}/idba'.format(threads,samplename,samplename))
        rename_contig = ('cp temp{}/idba/contig.fa {}.contigs.fasta'.format(samplename,samplename))
        fout.write('{}\n{}\n{}\n{}\n'.format(mk_temp_idba,merge_fq,assemble_fa,rename_contig))

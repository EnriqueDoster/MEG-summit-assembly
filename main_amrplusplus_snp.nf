#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

if (params.help ) {
    return help()
}
if( params.host_index ) {
    host_index = Channel.fromPath(params.host_index).toSortedList()
    //if( host_index.isEmpty() ) return index_error(host_index)
}
if( params.host ) {
    host = file(params.host)
    if( !host.exists() ) return host_error(host)
}
if( params.amr ) {
    amr = file(params.amr)
    if( !amr.exists() ) return amr_error(amr)
}
if( params.adapters ) {
    adapters = file(params.adapters)
    if( !adapters.exists() ) return adapter_error(adapters)
}
if( params.annotation ) {
    annotation = file(params.annotation)
    if( !annotation.exists() ) return annotation_error(annotation)
}
if( params.snp_annotation  ) {
    snp_annotation = file(params.snp_annotation)
}
if( params.hmm_analysis_script ) {
    hmm_analysis_script = file(params.hmm_analysis_script)
}
if( params.hmm_snp_annotation ) {
    hmm_snp_annotation = file(params.hmm_snp_annotation)
}
if( params.hmm_annotation ) {
    hmm_annotation = file(params.hmm_annotation)
}
if( params.hmm_group1 ) {
    hmm_group1 = file(params.hmm_group1)
}
if( params.hmm_group2 ) {
    hmm_group2 = file(params.hmm_group2)
}
if( params.hmm_group3 ) {
    hmm_group3 = file(params.hmm_group3)
}


threads = params.threads
smem_threads = params.smem_threads
threshold = params.threshold

min = params.min
max = params.max
skip = params.skip
samples = params.samples

leading = params.leading
trailing = params.trailing
slidingwindow = params.slidingwindow
minlen = params.minlen

Channel
    .fromFilePairs( params.reads, flat: true )
    .ifEmpty { exit 1, "Read pair files could not be found: ${params.reads}" }
    .set { reads }

process RunQC {
    tag { sample_id }


    publishDir "${params.output}/RunQC", mode: 'copy', pattern: '*.fastq',
        saveAs: { filename ->
            if(filename.indexOf("P.fastq") > 0) "Paired/$filename"
            else if(filename.indexOf("U.fastq") > 0) "Unpaired/$filename"
            else {}
        }

    input:
        set sample_id, file(forward), file(reverse) from reads

    output:
        set sample_id, file("${sample_id}.1P.fastq"), file("${sample_id}.2P.fastq") into (paired_fastq)
        set sample_id, file("${sample_id}.1U.fastq"), file("${sample_id}.2U.fastq") into (unpaired_fastq)
        file("${sample_id}.trimmomatic.stats.log") into (trimmomatic_stats)

    """
    ${JAVA} -jar ${TRIMMOMATIC}/trimmomatic-0.36.jar \
      PE \
      -threads ${threads} \
      $forward $reverse -baseout ${sample_id} \
      ILLUMINACLIP:${adapters}:2:30:10:3:TRUE \
      LEADING:${leading} \
      TRAILING:${trailing} \
      SLIDINGWINDOW:${slidingwindow} \
      MINLEN:${minlen} \
      2> ${sample_id}.trimmomatic.stats.log
    mv ${sample_id}_1P ${sample_id}.1P.fastq
    mv ${sample_id}_2P ${sample_id}.2P.fastq
    mv ${sample_id}_1U ${sample_id}.1U.fastq
    mv ${sample_id}_2U ${sample_id}.2U.fastq
    """
}

trimmomatic_stats.toSortedList().set { trim_stats }

process QCStats {
    tag { sample_id }

    publishDir "${params.output}/RunQC", mode: 'copy',
        saveAs: { filename ->
            if(filename.indexOf(".stats") > 0) "Stats/$filename"
            else {}
        }

    input:
        file(stats) from trim_stats

    output:
	file("trimmomatic.stats")

    """
    python3 $baseDir/bin/trimmomatic_stats.py -i ${stats} -o trimmomatic.stats
    """
}

if( !params.host_index ) {
    process BuildHostIndex {
        publishDir "${params.output}/BuildHostIndex", mode: "copy"

        tag { host.baseName }

        input:
            file(host)

        output:
            file '*' into (host_index)

        """
        bwa index ${host}
        """
    }
}

process AlignReadsToHost {
    tag { sample_id }

    publishDir "${params.output}/AlignReadsToHost", mode: "copy"

    input:
        set sample_id, file(forward), file(reverse) from paired_fastq
        file index from host_index.first()
        file host

    output:
        set sample_id, file("${sample_id}.host.sorted.bam") into (host_bam)

    """
    bwa mem ${host} ${forward} ${reverse} -t ${threads} > ${sample_id}.host.sam
    samtools view -bS ${sample_id}.host.sam | samtools sort -@ ${threads} -o ${sample_id}.host.sorted.bam
    rm *.sam
    """
}

process RemoveHostDNA {
    tag { sample_id }

    publishDir "${params.output}/RemoveHostDNA", mode: "copy", pattern: '*.bam',
	saveAs: { filename ->
            if(filename.indexOf(".bam") > 0) "NonHostBAM/$filename"
        }

    input:
        set sample_id, file(bam) from host_bam

    output:
        set sample_id, file("${sample_id}.host.sorted.removed.bam") into (non_host_bam)
        file("${sample_id}.samtools.idxstats") into (idxstats_logs)

    """
    samtools index ${bam} && samtools idxstats ${bam} > ${sample_id}.samtools.idxstats
    samtools view -h -f 4 -b ${bam} -o ${sample_id}.host.sorted.removed.bam
    """
}

idxstats_logs.toSortedList().set { host_removal_stats }

process HostRemovalStats {
    tag { sample_id }

    publishDir "${params.output}/RemoveHostDNA", mode: "copy",
        saveAs: { filename ->
            if(filename.indexOf(".stats") > 0) "HostRemovalStats/$filename"
        }

    input:
        file(stats) from host_removal_stats

    output:
        file("host.removal.stats")

    """
    python3 $baseDir/bin/samtools_idxstats.py -i ${stats} -o host.removal.stats
    """
}


if( !params.amr_index ) {
    process BuildAMRIndex {
        tag { amr.baseName }

        input:
            file(amr)

        output:
            file '*' into (amr_index)

        """
        bwa index ${amr}
        """
    }
}

process BAMToFASTQ {
    tag { sample_id }

    publishDir "${params.output}/BAMToFASTQ", mode: "copy"

    input:
        set sample_id, file(bam) from non_host_bam

    output:
        set sample_id, file("${sample_id}.non.host.R1.fastq.gz"), file("${sample_id}.non.host.R2.fastq.gz") into (non_host_fastq_dedup,non_host_fastq_assembly,non_host_fastq_megares,non_host_fastq_nodedup)

    """
    bedtools  \
       bamtofastq \
      -i ${bam} \
      -fq ${sample_id}.non.host.R1.fastq \
      -fq2 ${sample_id}.non.host.R2.fastq
    gzip *.fastq
    """
}

/* From here on you have the non-host reads for actual analysis */

/*
-
--
---
---- nonhost reads for megares
---
--
-
*/


process AlignToAMR {
     tag { sample_id }

     publishDir "${params.output}/AlignToAMR", mode: "copy"

     input:
         set sample_id, file(forward), file(reverse) from non_host_fastq_megares
         file index from amr_index.first()
         file amr

     output:
         set sample_id, file("${sample_id}.amr.alignment.sam") into (megares_resistome_sam, megares_rarefaction_sam, megares_snp_sam , megares_snpfinder_sam)
         set sample_id, file("${sample_id}.amr.alignment.dedup.sam") into (megares_dedup_resistome_sam)
         set sample_id, file("${sample_id}.amr.alignment.dedup.bam") into (megares_dedup_resistome_bam)
         

     """
     bwa mem ${amr} ${forward} ${reverse} -t ${threads} -R '@RG\\tID:${sample_id}\\tSM:${sample_id}' > ${sample_id}.amr.alignment.sam
     samtools view -S -b ${sample_id}.amr.alignment.sam > ${sample_id}.amr.alignment.bam
     samtools sort -n ${sample_id}.amr.alignment.bam -o ${sample_id}.amr.alignment.sorted.bam
     samtools fixmate ${sample_id}.amr.alignment.sorted.bam ${sample_id}.amr.alignment.sorted.fix.bam
     samtools sort ${sample_id}.amr.alignment.sorted.fix.bam -o ${sample_id}.amr.alignment.sorted.fix.sorted.bam
     samtools rmdup -S ${sample_id}.amr.alignment.sorted.fix.sorted.bam ${sample_id}.amr.alignment.dedup.bam
     samtools view -h -o ${sample_id}.amr.alignment.dedup.sam ${sample_id}.amr.alignment.dedup.bam
     rm ${sample_id}.amr.alignment.bam
     rm ${sample_id}.amr.alignment.sorted*.bam
     """
}

process RunResistome {
    tag { sample_id }

    publishDir "${params.output}/RunResistome", mode: "copy"

    input:
        set sample_id, file(sam) from megares_resistome_sam
        file annotation
        file amr

    output:
        file("${sample_id}.gene.tsv") into (megares_resistome_counts, SNP_confirm_long)

    """
    resistome -ref_fp ${amr} \
      -annot_fp ${annotation} \
      -sam_fp ${sam} \
      -gene_fp ${sample_id}.gene.tsv \
      -group_fp ${sample_id}.group.tsv \
      -class_fp ${sample_id}.class.tsv \
      -mech_fp ${sample_id}.mechanism.tsv \
      -t ${threshold}
    """
}

megares_resistome_counts.toSortedList().set { megares_amr_l_to_w }

process AMRLongToWide {
    tag { }

    publishDir "${params.output}/AMRLongToWide", mode: "copy"

    input:
        file(resistomes) from megares_amr_l_to_w

    output:
        file("AMR_analytic_matrix.csv") into amr_master_matrix

    """
    mkdir ret
    python3 $baseDir/bin/amr_long_to_wide.py -i ${resistomes} -o ret
    mv ret/AMR_analytic_matrix.csv .
    """
}


/* samtools deduplication of megares alignment */
process SamDedupRunResistome {
    tag { sample_id }

    publishDir "${params.output}/SamDedupRunResistome", mode: "copy"

    input:
        set sample_id, file(sam) from megares_dedup_resistome_sam
        file annotation
        file amr

    output:
        file("${sample_id}.gene.tsv") into (megares_dedup_resistome_counts)

    """
    resistome -ref_fp ${amr} \
      -annot_fp ${annotation} \
      -sam_fp ${sam} \
      -gene_fp ${sample_id}.gene.tsv \
      -group_fp ${sample_id}.group.tsv \
      -class_fp ${sample_id}.class.tsv \
      -mech_fp ${sample_id}.mechanism.tsv \
      -t ${threshold}
    """
}

megares_dedup_resistome_counts.toSortedList().set { megares_dedup_amr_l_to_w }

process SamDedupAMRLongToWide {
    tag { }

    publishDir "${params.output}/SamDedup_AMRLongToWide", mode: "copy"

    input:
        file(resistomes) from megares_dedup_amr_l_to_w

    output:
        file("SamDedup_AMR_analytic_matrix.csv") into megares_dedup_amr_master_matrix

    """
    mkdir ret
    python3 $baseDir/bin/amr_long_to_wide.py -i ${resistomes} -o ret
    mv ret/AMR_analytic_matrix.csv SamDedup_AMR_analytic_matrix.csv
    """
}




/*
-
--
---
---- SNP analysis
---
--
-
*/




process RunFreebayes {
    tag { sample_id }

    publishDir "${params.output}/RunFreebayes", mode: "copy"

    input:
        set sample_id, file(bam) from megares_dedup_resistome_bam
        file annotation
        file amr

    output:
        set sample_id, file("${sample_id}.results.vcf") into (SNP)

    """
    freebayes \
      -f ${amr} \
      -p 1 \
      ${bam} > ${sample_id}.results.vcf
    """
}


process RunRarefaction {
    tag { sample_id }

    publishDir "${params.output}/RunRarefaction", mode: "copy"

    input:
        set sample_id, file(sam) from megares_rarefaction_sam
        file annotation
        file amr

    output:
        set sample_id, file("*.tsv") into (rarefaction)

    """
    rarefaction \
      -ref_fp ${amr} \
      -sam_fp ${sam} \
      -annot_fp ${annotation} \
      -gene_fp ${sample_id}.gene.tsv \
      -group_fp ${sample_id}.group.tsv \
      -class_fp ${sample_id}.class.tsv \
      -mech_fp ${sample_id}.mech.tsv \
      -min ${min} \
      -max ${max} \
      -skip ${skip} \
      -samples ${samples} \
      -t ${threshold}
    """
}

process RunSNPFinder {
    tag { sample_id }

    publishDir "${params.output}/RunSNPFinder", mode: "copy"

    input:
        set sample_id, file(sam) from megares_snpfinder_sam
        file amr

    output:
        set sample_id, file("*.tsv") into (snp)

    """
    snpfinder \
      -amr_fp ${amr} \
      -sampe ${sam} \
      -out_fp ${sample_id}.tsv
    """
}



def nextflow_version_error() {
    println ""
    println "This workflow requires Nextflow version 0.25 or greater -- You are running version $nextflow.version"
    println "Run ./nextflow self-update to update Nextflow to the latest available version."
    println ""
    return 1
}

def adapter_error(def input) {
    println ""
    println "[params.adapters] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}

def amr_error(def input) {
    println ""
    println "[params.amr] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}

def annotation_error(def input) {
    println ""
    println "[params.annotation] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}

def fastq_error(def input) {
    println ""
    println "[params.reads] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}

def host_error(def input) {
    println ""
    println "[params.host] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}

def index_error(def input) {
    println ""
    println "[params.host_index] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}

def help() {
    println ""
    println "Program: AmrPlusPlus"
    println "Documentation: https://github.com/colostatemeg/amrplusplus/blob/master/README.md"
    println "Contact: Christopher Dean <cdean11@colostate.edu>"
    println ""
    println "Usage:    nextflow run main.nf [options]"
    println ""
    println "Input/output options:"
    println ""
    println "    --reads         STR      path to FASTQ formatted input sequences"
    println "    --adapters      STR      path to FASTA formatted adapter sequences"
    println "    --host          STR      path to FASTA formatted host genome"
    println "    --host_index    STR      path to BWA generated index files"
    println "    --amr           STR      path to AMR resistance database"
    println "    --annotation    STR      path to AMR annotation file"
    println "    --output        STR      directory to write process outputs to"
    println "    --KRAKENDB      STR      path to kraken database"
    println ""
    println "Trimming options:"
    println ""
    println "    --leading       INT      cut bases off the start of a read, if below a threshold quality"
    println "    --minlen        INT      drop the read if it is below a specified length"
    println "    --slidingwindow INT      perform sw trimming, cutting once the average quality within the window falls below a threshold"
    println "    --trailing      INT      cut bases off the end of a read, if below a threshold quality"
    println ""
    println "Algorithm options:"
    println ""
    println "    --threads       INT      number of threads to use for each process"
    println "    --threshold     INT      gene fraction threshold"
    println "    --min           INT      starting sample level"
    println "    --max           INT      ending sample level"
    println "    --samples       INT      number of sampling iterations to perform"
    println "    --skip          INT      number of levels to skip"
    println ""
    println "Help options:"
    println ""
    println "    --help                   display this message"
    println ""
    return 1
}

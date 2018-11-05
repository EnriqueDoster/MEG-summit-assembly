#!/usr/bin/env nextflow

threads = params.threads




if (params.help ) {
    return help()
}
if( !nextflow.version.matches('0.25+') ) {
    return nextflow_version_error()
}
if( params.host_index ) {
    host_index = Channel.fromPath(params.host_index).toSortedList()
    //if( host_index.isEmpty() ) return index_error(host_index)
}
if( params.host ) {
    host = file(params.host)
    if( !host.exists() ) return host_error(host)
}
if( params.adapters ) {
    adapters = file(params.adapters)
    if( !adapters.exists() ) return adapter_error(adapters)
}



// Trimmomatic options
leading = params.leading
trailing = params.trailing
slidingwindow = params.slidingwindow
minlen = params.minlen
adapters = params.adapters


Channel
    .fromFilePairs( params.reads, flat: true )
    .ifEmpty { exit 1, "Read pair files could not be found: ${params.reads}" }
    .into { read_pairs; fastqc_pairs }

process QualityControl {
    tag { dataset_id }

    publishDir "${params.output}/QualityControlOutput", mode: "symlink",
        saveAs: { filename ->
            if(filename.indexOf("P.fastq") > 0) "Paired/$filename"
            else if(filename.indexOf(".log") > 0) "Log/$filename"
            else {}
        }

    input:
        set dataset_id, file(forward), file(reverse) from read_pairs

    output:
        set dataset_id, file("${dataset_id}.1P.fastq.gz"), file("${dataset_id}.2P.fastq.gz") into (paired_fastq)
        set dataset_id, file("${dataset_id}.trimmomatic.stats.log") into (trimmomatic_logs)

    """
    /usr/lib/jvm/java-7-openjdk-amd64/bin/java -jar ${TRIMMOMATIC}/trimmomatic-0.36.jar \
        PE \
        -threads ${threads} \
        $forward $reverse -baseout ${dataset_id} \
        ILLUMINACLIP:${adapters}:2:30:10:3:TRUE \
        LEADING:${leading} \
        TRAILING:${trailing} \
        SLIDINGWINDOW:${slidingwindow} \
        MINLEN:${minlen} \
        2> ${dataset_id}.trimmomatic.stats.log

    gzip -c ${dataset_id}_1P > ${dataset_id}.1P.fastq.gz
    gzip -c ${dataset_id}_2P > ${dataset_id}.2P.fastq.gz
    rm ${dataset_id}_1P
    rm ${dataset_id}_2P
    rm ${dataset_id}_1U
    rm ${dataset_id}_2U
    """
}

/*
trimmomatic_stats.toSortedList().set { trim_stats }

process QCStats {
    tag { dataset_id }

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
*/


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
    tag { dataset_id }

    publishDir "${params.output}/AlignReadsToHost", mode: "copy"

    input:
        set dataset_id, file(forward), file(reverse) from paired_fastq
        file index from host_index.first()
        file host

    output:
        set dataset_id, file("${dataset_id}.host.sam") into (host_sam)

    """
    bwa mem ${host} ${forward} ${reverse} -t ${threads} > ${dataset_id}.host.sam
    """
}

process RemoveHostDNA {
    tag { dataset_id }

    publishDir "${params.output}/RemoveHostDNA", mode: "copy", pattern: '*.bam',
	saveAs: { filename ->
            if(filename.indexOf(".bam") > 0) "NonHostBAM/$filename"
        }

    input:
        set dataset_id, file(sam) from host_sam

    output:
        set dataset_id, file("${dataset_id}.host.sorted.removed.bam") into (non_host_bam)
        file("${dataset_id}.samtools.idxstats") into (idxstats_logs)

    """
    samtools view -bS ${sam} | samtools sort -@ ${threads} -o ${dataset_id}.host.sorted.bam
    samtools index ${dataset_id}.host.sorted.bam && samtools idxstats ${dataset_id}.host.sorted.bam > ${dataset_id}.samtools.idxstats
    samtools view -h -f 4 -b ${dataset_id}.host.sorted.bam -o ${dataset_id}.host.sorted.removed.bam
    """
}

idxstats_logs.toSortedList().set { host_removal_stats }

/*
process HostRemovalStats {
    tag { dataset_id }

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
*/
process BAMToFASTQ {
    tag { dataset_id }

    publishDir "${params.output}/BAMToFASTQ", mode: "copy"

    input:
        set dataset_id, file(bam) from non_host_bam

    output:
        set dataset_id, file("${dataset_id}.non.host.R1.fastq"), file("${dataset_id}.non.host.R2.fastq") into (non_host_fastq, non_host_fastq_kraken)

    """
    bedtools  \
       bamtofastq \
      -i ${bam} \
      -fq ${dataset_id}.non.host.R1.fastq \
      -fq2 ${dataset_id}.non.host.R2.fastq
    """
}
/*
process AssembleReads {
    tag { dataset_id }

    publishDir "${params.output}/AssembledFiles", mode: "symlink"

    input:
        set dataset_id, file(forward), file(reverse) from non_host_fastq

    output:
        set dataset_id, file("${dataset_id}.contigs.fa") into (idba_assemblies)

    script:
    """
    mkdir -p temp/idba
    fq2fa --merge --filter <( zcat $forward) <( zcat $reverse ) temp/interleavened.fasta
    idba_ud --num_threads 60 -l temp/interleavened.fasta -o temp/idba

	cp temp/idba/contig.fa ${dataset_id}.contigs.fasta

    """
}
*/

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
    println "Program: summit-assembly"
    println "Version: $workflow.repository - $workflow.revision [$workflow.commitId]"
    println "Contact: Steven Lakin <steven.m.lakin@gmail.com>"
    println ""
    println "Usage:    nextflow run summit-assembly.nf [options]"
    println ""
    println "Input/output options:"
    println ""
    println "    --reads         STR      path to FASTQ formatted input sequences"
    println "    --output        STR      directory to write process outputs to"
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
    println ""
    println "Help options:"
    println ""
    println "    --help                   display this message"
    println ""
    return 1
}

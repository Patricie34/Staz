nextflow.enable.dsl=2


process test {
  container = "staphb/fastqc:latest" 

  script:
  """
sleep infinity 
  """
}

// Run FASTQC
process fastqc {
publishDir("${params.outdir}/fastqc", mode: 'copy')

  input:
  tuple val(sample_id), path(reads)

  output:
  path "fastqc_${sample_id}/*"

  container "staphb/fastqc:latest" 

  script:
  """
  mkdir fastqc_${sample_id}
  fastqc -o fastqc_${sample_id}  ${reads}
  """
}


// Run BWA
process bwa {
    //conda '/opt/miniconda/envs/projekt'
    memory "7 GB"
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}.sam")

    container "biocontainers/bwa:v0.7.17_cv1"

    script:
    """
    bwa mem -R "@RG\\tID:HG002_HiSeq\\tSM:HG002_HiSeq\\tPL:illumina" $params.genome \
    $reads > ${sample_id}.sam
    """
}

// Convert .sam file to .bam
process sam_to_bam {
    publishDir("${params.outdir}/${sample_id}/bam", mode: 'copy')
    memory '20 MB'

    input:
    tuple val(sample_id), path(sam)

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    container "staphb/samtools:latest"

    script:
    """
    samtools view -bS $sam \
    > ${sample_id}.bam
    """
}

// Sort .bam files
process sort_bam {
    publishDir("${params.outdir}/${sample_id}/sort_bam", mode: 'copy')
    memory '50 MB'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam")

    container "staphb/samtools:latest"

    script:
    """
    samtools sort $bam \
    > ${sample_id}_sorted.bam
    """
}

// Extract only mapped reads
process extract_mapped_reads {
    publishDir("${params.outdir}/${sample_id}/sort_mapped_bam", mode: 'copy')
    memory '20 MB'

    input:
    tuple val(sample_id), path(sort_bam)

    output:
    tuple val(sample_id), path("${sample_id}_sorted_mapped.bam")

    container "staphb/samtools:latest"

    script:
    """
   samtools view -bS -F 4 $sort_bam \
    > ${sample_id}_sorted_mapped.bam
    """
}

// // Indexing the sorted and mapped .bam file
// process index_bam {
//     publishDir("${params.outdir}/${sample_id}/indexed_bam", mode: 'copy')
//     input:
//     tuple val(sample_id), path(sort_mapped_bam)

//     output:
//     // tuple val(sample_id), path(sort_mapped_bam), path("${sort_mapped_bam}.bai") //zpusob 1
//     tuple val(sample_id), path("${sort_mapped_bam}.bai") //zpusob 2 (musi byt nastaven Channel BaiBam)

//     container "staphb/samtools:latest"

//     script:
//     """
//     samtools index $sort_mapped_bam
//     """
// }

// Variant calling

// Mark Duplicates
process mark_duplicates {
    memory '10 GB'
    publishDir("${params.outdir}/${sample_id}/MarkDuplicates", mode: 'copy')
    input:
    tuple val(sample_id), path(bam) //path(bai), path(bam)
    
    output:
    tuple val(sample_id), path("${sample_id}_sorted_mapped_MD.bam"), path("${sample_id}_marked_dup_metrics.txt")

    container "mgibio/gatk-cwl:4.1.8.1"

    script:
    """
    java -jar $params.gatk_package MarkDuplicates \
        -I "$bam" \
        -O "${sample_id}_sorted_mapped_MD.bam" \
        -M "${sample_id}_marked_dup_metrics.txt"
    """
}

// Indexing the sorted and marked duplicates BAM file
process index_MD_bam {
    publishDir("${params.outdir}/${sample_id}/indexed_bam", mode: 'copy')
    memory '50 MB'

    input:
    tuple val(sample_id), path(sorted_mapped_MD_bam)

    output:
    // tuple val(sample_id), path(sort_mapped_bam), path("${sort_mapped_bam}.bai") //zpusob 1
    tuple val(sample_id), path("${sorted_mapped_MD_bam}.bai") //zpusob 2 (musi byt nastaven Channel BaiBam)

    container "staphb/samtools:latest"

    script:
    """
    samtools index $sorted_mapped_MD_bam
    """
}

// Base Recalibrator
process base_recal {
    memory '1 GB'

    input:
    tuple val(sample_id), path(sorted_mapped_MD_bai), path(sorted_mapped_MD_bam)

    output:
    tuple val(sample_id), path("${sample_id}_sorted_mapped_recal_data.table")
    
    container "mgibio/gatk-cwl:4.1.8.1"
    
    script:
    """
    java -jar $params.gatk_package BaseRecalibrator \
   -I $sorted_mapped_MD_bam \
   -R $params.ref_genome \
   --known-sites $params.known_sites \
   -O ${sample_id}_sorted_mapped_recal_data.table
   """
}

// ApplyBQSR
process apply_bqsr {
    memory '1 GB'
    input:
    tuple val(sample_id), path(sorted_mapped_recal_data_table), path(sort_MD_bai), path(sort_MD_bam)

    output:
    tuple val(sample_id), path("${sample_id}_sorted_mapped_recal.bam")

    container "mgibio/gatk-cwl:4.1.8.1"

    script:
    """
    java -jar $params.gatk_package ApplyBQSR \
   -R $params.ref_genome \
   -I $sort_MD_bam \
   --bqsr-recal-file $sorted_mapped_recal_data_table \
   -O ${sample_id}_sorted_mapped_recal.bam
   """
}

// Index Recalibrate BAMs
process index_recal_bam {
	memory '20 MB'

	input:
	tuple val(sample_id), path(sorted_mapped_recal_bam)

	output:
	tuple val(sample_id), path("${sample_id}_sorted_mapped_recal.bam"), path("${sample_id}_sorted_mapped_recal.bam.bai")

	script:
	"""
	samtools index $sorted_mapped_recal_bam
	"""
}
// samtools merge --threads ${sample_id}_sorted_mapped_recal.bam ${bams}

// HaplotypeCaller
process hapl_caller {
    memory "3 GB"
    publishDir("${params.outdir}/${sample_id}/hapl_caller", mode: 'copy')

    input:
    tuple val(sample_id), path(sorted_mapped_recal_bam), path(sorted_mapped_recal_bam_bai)

    output:
    tuple val(sample_id), path("${sample_id}_sorted_mapped_recal_haplCall.g.vcf.gz")

    container "mgibio/gatk-cwl:4.1.8.1"

    script:
    """
    java -jar $params.gatk_package HaplotypeCaller \
   -R $params.ref_genome \
   -I $sorted_mapped_recal_bam \
   -O ${sample_id}_sorted_mapped_recal_haplCall.g.vcf.gz \
   -ERC GVCF
   """
}

// -Xmx4g

// Benchmark
process benchmark {
    memory "100 GB"
    publishDir("${params.outdir}/benchmark", mode: 'copy')

    input:
    tuple val(sample_id), path(mgi), path(ilumina), path(true_vcf)
    // path(true_vcf), path(mgi_vcf), path(ilumina_vcf) //zpusob 1
    
    output:
    path("*")

    container "biocontainers/biopython:v1.73dfsg-1-deb-py3_cv1"

    script:
    """
    python $params.script --true ${true_vcf} --mgi ${mgi} --ilumina ${ilumina}

   """
   //python $params.script --true ${true_vcf} --mgi ${mgi} --ilumina ${ilumina}
    //  python benchmark.py --true ${true_vcf} --mgi ${mgi_vcf} --ilumina ${ilumina_vcf} // zpusob 1
    //python benchmark.py --true ${params.true_vcf} --mgi ${params.mgi_vcf} --ilumina ${params.ilumina_vcf}  zpusob 2
// echo TEST
//      python $params.script --true ${params.true_vcf} --mgi ${params.mgi_vcf} --ilumina ${params.ilumina_vcf} 
}



// Define workflow
workflow {
    //  Channels
// ref_ch = Channel.fromPath(params.genome, checkIfExists: true)
// ref2_ch = Channel.fromPath(params.ref_genome, checkIfExists: true)
reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)//.view()
true_vcf_ch = Channel.fromPath(params.true_vcf, checkIfExists: true)//.view()
//test()
//reads_ch.view()
// fastqc(reads_ch)                                                                                                                                                                                                                                                                    
sam_ch = bwa(reads_ch)
bam_ch = sam_to_bam(sam_ch)
sort_ch = sort_bam(bam_ch)
mapp_reads = extract_mapped_reads(sort_ch)
mark_ch = mark_duplicates(mapp_reads)//.view()
mark_bam_ch = mark_ch.map{it -> [it[0], it[1]]}//.view()
indexed_mark_ch = index_MD_bam(mark_bam_ch)//.view()
BaiBam_ch = indexed_mark_ch.join(mark_bam_ch)//.view() //zpusob 2
base_ch = base_recal(BaiBam_ch)//.view()
MD_base_ch = base_ch.join(BaiBam_ch)//.view()
recal_ch = apply_bqsr(MD_base_ch)//.view()
recal_index_ch = index_recal_bam(recal_ch)//.view()
hapl_ch = hapl_caller(recal_index_ch)
hapl_caller_ch = hapl_ch.map{it -> [it[0][0..5], it[1]]}.groupTuple()//.view()
call_ch = hapl_caller_ch.combine(true_vcf_ch).map{it -> [it[0], it[1][0], it[1][1], it[2]]}//.view()
benchmark(call_ch)
}
    

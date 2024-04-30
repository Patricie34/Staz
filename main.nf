nextflow.enable.dsl=2

// Run FASTQC
process fastqc {
publishDir("${params.outdir}/fastqc", mode: 'copy')

  input:
  tuple val(sample_id), path(reads)

  output:
  path "fastqc_${sample_id}/*"

  container = "staphb/fastqc:latest" 

  script:
  """
sleep infinity 
  mkdir fastqc_${sample_id}
  fastqc -o fastqc_${sample_id}  ${reads}
  """
}


// Run BWA
process bwa {
    conda '/opt/miniconda/envs/projekt'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}.sam")

    script:
    """
    bwa mem -R "@RG\\tID:HG002_HiSeq\\tSM:HG002_HiSeq\\tPL:illumina" $params.genome \
    $reads > ${sample_id}.sam
    """
}

// Convert .sam file to .bam
process sam_to_bam {
    publishDir("${params.outdir}/${sample_id}/bam", mode: 'copy')
    memory '100 MB'

    input:
    tuple val(sample_id), path(sam)

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    script:
    """
    samtools view -bS $sam \
    > ${sample_id}.bam
    """
}

// Sort .bam files
process sort_bam {
    publishDir("${params.outdir}/${sample_id}/sort_bam", mode: 'copy')
    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam")

    script:
    """
    samtools sort $bam \
    > ${sample_id}_sorted.bam
    """
}

// Extract only mapped reads
process extract_mapped_reads {
    publishDir("${params.outdir}/${sample_id}/sort_mapped_bam", mode: 'copy')
    input:
    tuple val(sample_id), path(sort_bam)

    output:
    tuple val(sample_id), path("${sample_id}_sorted_mapped.bam")

    script:
    """
   samtools view -bS -F 4 $sort_bam \
    > ${sample_id}_sorted_mapped.bam
    """
}

// Indexing the sorted and mapped .bam file
process index_bam {
    publishDir("${params.outdir}/${sample_id}/indexed_bam", mode: 'copy')
    input:
    tuple val(sample_id), path(sort_mapped_bam)

    output:
    // tuple val(sample_id), path(sort_mapped_bam), path("${sort_mapped_bam}.bai") //zpusob 1
    tuple val(sample_id), path("${sort_mapped_bam}.bai") //zpusob 2 (musi byt nastaven Channel BaiBam)

    script:
    """
    samtools index $sort_mapped_bam
    """
}

// Variant calling

// Mark Duplicates
process mark_duplicates {
    publishDir("${params.outdir}/${sample_id}/MarkDuplicates", mode: 'copy')
    input:
    tuple val(sample_id), path(bai), path(bam)
    
    output:
    tuple val(sample_id), path("${sample_id}_sorted_mapped_MD.bam"), path("${sample_id}_marked_dup_metrics.txt")

    script:
    """
    java -jar $params.gatk_package MarkDuplicates \
        -I "$bam" \
        -O "${sample_id}_sorted_mapped_MD.bam" \
        -M "${sample_id}_marked_dup_metrics.txt"
    """
}

// Base Recalibrator
process base_recal {
    input:
    tuple val(sample_id), path(sorted_mapped_MD_bam)

    output:
    tuple val(sample_id), path("${sample_id}_sorted_mapped_recal_data.table")

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
    input:
    tuple val(sample_id), path(sorted_mapped_recal_data_table), path(sort_MD_bam)

    output:
    tuple val(sample_id), path("${sample_id}_sorted_mapped_recal.bam")

    script:
    """
    java -jar $params.gatk_package ApplyBQSR \
   -R $params.ref_genome \
   -I $sort_MD_bam \
   --bqsr-recal-file $sorted_mapped_recal_data_table \
   -O ${sample_id}_sorted_mapped_recal.bam
   """
}

// HaplotypeCaller
process hapl_caller {
    publishDir("${params.outdir}/${sample_id}/hapl_caller", mode: 'copy')

    input:
    tuple val(sample_id), path(sorted_mapped_recal_bam)

    output:
    tuple val(sample_id), path("${sample_id}_sorted_mapped_recal_haplCall.g.vcf.gz")

    script:
    """
    java -Xmx4g -jar $params.gatk_package HaplotypeCaller \
   -R $params.ref_genome \
   -I $sorted_mapped_recal_bam \
   -O ${sample_id}_sorted_mapped_recal_haplCall.g.vcf.gz \
   -ERC GVCF
   """
}

// Benchmark
process benchmark {
    publishDir("${params.outdir}/benchmark", mode: 'copy')

    input:
    tuple val(sample_id), path(mgi), path(ilumina), path(true_vcf)
    // path(true_vcf), path(mgi_vcf), path(ilumina_vcf) //zpusob 1
    
    output:
    path("*")

    script:
    """
    python $params.script --true ${true_vcf} --mgi ${mgi} --ilumina ${ilumina}
     
   """
    //  python benchmark.py --true ${true_vcf} --mgi ${mgi_vcf} --ilumina ${ilumina_vcf} // zpusob 1
    //python benchmark.py --true ${params.true_vcf} --mgi ${params.mgi_vcf} --ilumina ${params.ilumina_vcf}  zpusob 2
// echo TEST
//      python $params.script --true ${params.true_vcf} --mgi ${params.mgi_vcf} --ilumina ${params.ilumina_vcf} 
}



// Define workflow
workflow {
    //  Channels
//ref_ch = Channel.fromPath(params.genome, checkIfExists: true)
//ref2_ch = Channel.fromPath(params.ref_genome, checkIfExists: true)
reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)//.view()
true_vcf_ch = Channel.fromPath(params.true_vcf, checkIfExists: true)//.view()

//reads_ch.view()
    fastqc(reads_ch)
//     sam_ch = bwa(reads_ch)
// //sam_ch.view()
//     bam_ch = sam_to_bam(sam_ch)
//     sort_ch = sort_bam(bam_ch)
//     mapp_reads = extract_mapped_reads(sort_ch)
//     indexed_ch = index_bam(mapp_reads)//.view()
//     //indexed_ch.view()
//     BaiBam_ch = indexed_ch.join(mapp_reads)//.view() //zpusob 2
//     //mark_duplicates(BaiBam_ch).view()
//     mark_ch = mark_duplicates(BaiBam_ch)//.view()
//     mark_bam_ch = mark_ch.map{it -> [it[0], it[1]]}//.view()
//     base_ch = base_recal(mark_bam_ch)
//     MD_base_ch = base_ch.join(mark_bam_ch)//.view()
//     //apply_bqsr(MD_base_ch)//.view()
//     recal_ch = apply_bqsr(MD_base_ch)
//     //hapl_caller(recal_ch)//.view()
//     //mark_ch.map{it -> [it[0][0..5], it[1]]}.view()
//     hapl_ch = hapl_caller(recal_ch)
//     hapl_caller_ch = hapl_ch.map{it -> [it[0][0..5], it[1]]}.groupTuple()//.view()
//     call_ch = hapl_caller_ch.combine(true_vcf_ch).map{it -> [it[0], it[1][0], it[1][1], it[2]]}.view()
//     benchmark(call_ch)
}
    

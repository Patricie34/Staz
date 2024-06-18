nextflow.enable.dsl=2

process bcftools{

    memory '20 MB'

    input:
    tuple val(sample_id), path(bcf), path(bcf_sci)

    output:
    tuple val(sample_id), path("${sample_id}.sv.tsv"), path("${sample_id}.input.tsv")

    container "staphb/bcftools:latest"

    script:
    """
	bcftools query -f "%CHROM\t%POS\t%CHROM\t%INFO/END\t%ID\n" ${bcf} | grep -v "BND" > ${sample_id}.sv.tsv 
	bcftools query -f "%CHROM\t%POS\t%INFO/CHR2\t%INFO/POS2\t%ID\n" ${bcf} | grep "BND" >> ${sample_id}.sv.tsv
    cat ${sample_id}.sv.tsv | cut -f 1,2,5 | sed 's/^chr//' | awk '{print \$1"\\t"(\$2-100)"\\t"(\$2+100)"\\t"\$3"Left";}' > ${sample_id}.input.tsv
	cat ${sample_id}.sv.tsv | cut -f 3,4,5 | sed 's/^chr//' | awk '{print \$1"\\t"(\$2-100)"\\t"(\$2+100)"\\t"\$3"Right";}' >> ${sample_id}.input.tsv   
    """
}


process alfred{
    publishDir("/storage/01.NanoBreak/data/samples/${sample_id}/delly_hg38/", mode: 'copy')
    memory '70 MB'

    input:
    tuple val(sample_id), path("${sample_id}.sv.tsv"), path("${sample_id}.input.tsv")

    output:
    tuple val(sample_id), path("${sample_id}.sv.gene.tsv")

    container "trausch/alfred:latest"

    script:
    """
    alfred annotate -d 3000 -g /home/user/delly_k8s/Homo_sapiens.GRCh38.107.gtf.gz -o ${sample_id}.sv.gene.tsv ${sample_id}.input.tsv
	rm ${sample_id}.sv.tsv ${sample_id}.input.tsv gene.bed
    """
}

workflow {
    bcf_ch = Channel.fromFilePairs(["/home/user/delly_k8s/delly_allSamples/*_delly.bcf", "/home/user/delly_k8s/delly_allSamples/*_delly.bcf.csi"])
    | map { id, files ->
    (sample, replicate, type) = id.tokenize("_")
    meta = [sample:sample, replicate:replicate, type:type]
    [sample, files]
}
    | map{it -> [it[0], it[1][0], it[1][1]]}

    | view
    tsv_ch = bcftools(bcf_ch)
    alf_ch = alfred(tsv_ch) //(tsv_ch.first())

}


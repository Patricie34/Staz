nextflow.enable.dsl=2

process bcftools{

    publishDir("${params.outdir}/${sample_id}/bcftools", mode: 'copy')

    input:
    tuple val(sample_id), path(bcf), path(bcf_sci)

    output:
    tuple val(sample_id), path("${sample_id}.sv.tsv"), path("${sample_id}.input.tsv")

    container "staphb/bcftools:latest"

    script:
    """
	bcftools query -f "%CHROM\t%POS\t%CHROM\t%INFO/END\t%ID\n" ${bcf} | grep -v "BND" > ${sample_id}.sv.tsv 
	bcftools query -f "%CHROM\t%POS\t%INFO/CHR2\t%INFO/POS2\t%ID\n" ${bcf} | grep "BND" >> ${sample_id}}.sv.tsv
	cat ${sample_id}}.sv.tsv | cut -f 1,2,5 | sed 's/^chr//' | awk '{print $1"\t"($2-100)"\t"($2+100)"\t"$3"Left";}' > ${sample_id}}.input.tsv
	cat ${sample_id}}.sv.tsv | cut -f 3,4,5 | sed 's/^chr//' | awk '{print $1"\t"($2-100)"\t"($2+100)"\t"$3"Right";}' >> ${sample_id}}.input.tsv

    """
}

process alfred{
    input:

    output:

    script:
    """
    /opt/dev/alfred/src/alfred annotate -d 3000 -g /opt/dev/alfred/gtf/Homo_sapiens.GRCh38.107.gtf.gz -o ${tumorDNAIllumina}.sv.gene.tsv ${tumorDNAIllumina}.input.tsv
	rm ${tumorDNAIllumina}.sv.tsv ${tumorDNAIllumina}.input.tsv gene.bed
    """
}

// workflow {
//     Channel.fromFilePairs(["/home/user/delly_k8s/delly_allSamples/*_delly.bcf", "/home/user/delly_k8s/delly_allSamples/*_delly.bcf.csi"], checkIfExists: true)
//     | view
// }


// workflow {
//     Channel.fromFilePairs(["/home/user/delly_k8s/delly_allSamples/*_delly.bcf", "/home/user/delly_k8s/delly_allSamples/*_delly.bcf.csi"])
//     | map { id, reads ->
//         tokens = id.tokenize("_")
//     }
//     | view
// }


workflow {
    bcf_ch = Channel.fromFilePairs(["/home/user/delly_k8s/delly_allSamples/*_delly.bcf", "/home/user/delly_k8s/delly_allSamples/*_delly.bcf.csi"])
    | map { id, files ->
    (sample, replicate, type) = id.tokenize("_")
    meta = [sample:sample, replicate:replicate, type:type]
    [sample, files]
}
    tsv_ch = bcftools(bcf_ch)

}


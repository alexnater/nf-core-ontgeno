process BCFTOOLS_MPILEUP {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/474a5ea8dc03366b04df884d89aeacc4f8e6d1ad92266888e7a8e7958d07cde8/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f'}"

    input:
    tuple val(meta), path(bam), path(bai), path(intervals)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*vcf.gz"), path("*vcf.gz.tbi"), emit: vcf_tbi
    tuple val(meta), path("*stats.txt"), emit: stats
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def intervals_cmd = intervals ? "-R ${intervals}" : ""
    """
    bcftools \\
        mpileup \\
        --threads ${task.cpus} \\
        --fasta-ref ${fasta} \\
        ${intervals_cmd} \\
        ${args} \\
        ${bam} \\
        | bcftools call --threads ${task.cpus} ${args2} --output-type z --write-index=tbi --output ${prefix}.vcf.gz

    bcftools stats ${prefix}.vcf.gz > ${prefix}.bcftools_stats.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bcftools_stats.txt
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}

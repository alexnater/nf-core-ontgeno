process WHATSHAP_STATS {
    tag "${meta.id} - ${meta.caller}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/02/02705a968727688fe6e2c2943cb0967726f7cf32521c8f92d9349fd8509798c1/data':
        'community.wave.seqera.io/library/bcftools_htslib_whatshap:259a60305e8ff1a9' }"

    input:
    tuple val(meta), path(vcf), path(tbi), val(samples)

    output:
    tuple val(meta), path("*.gtf"), emit: gtf
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    for sample in ${samples.join(' ')}; do
        whatshap \\
            stats \\
            $args \\
            --sample \$sample \\
            --gtf ${prefix}.\${sample}.gtf \\
            $vcf
    
        sort -k1,1 -k4,4n ${prefix}.\${sample}.gtf | \\
            awk -v OFS="\\t" -v sample=\${sample} \\
                'BEGIN{counter=1} {print \$1,\$4-1,\$5,sample"."counter; ++counter}' \\
            >> ${prefix}.blocks.bed
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: \$(whatshap --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    for sample in ${samples.join(' ')}; do
        touch ${prefix}.\${sample}.gtf
    done

    touch ${prefix}.blocks.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: \$(whatshap --version)
    END_VERSIONS
    """
}
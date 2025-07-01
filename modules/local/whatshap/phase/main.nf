process WHATSHAP_PHASE {
    tag "${meta.id} - ${meta.caller}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/02/02705a968727688fe6e2c2943cb0967726f7cf32521c8f92d9349fd8509798c1/data':
        'community.wave.seqera.io/library/bcftools_htslib_whatshap:259a60305e8ff1a9' }"

    input:
    tuple val(meta),  path(variants)
    tuple val(meta2), path(reads), path(indices)
    tuple val(meta3), path(fasta), path(fai)

    output:
    tuple val(meta), path("${prefix}.phased.vcf.gz"), path("${prefix}.phased.vcf.gz.tbi"), emit: vcf_tbi
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args ?: ''
    def args3 = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def convert_bcf = variants.extension == 'bcf' ? 1 : 0

    """
    if [ $convert_bcf -eq 1 ]
    then
        bcftools \\
            view \\
            $args \\
            -o ${prefix}.vcf \\
            $variants
    else
        ln -s $variants ${prefix}.vcf
    fi

    whatshap \\
        phase \\
        $args2 \\
        --reference $fasta \\
        -o ${prefix}.phased.vcf.gz \\
        ${prefix}.vcf \\
        $reads

    tabix \\
        $args3 \\
        -p vcf \\
        ${prefix}.phased.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        whatshap: \$(whatshap --version)
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf
    touch ${prefix}.phased.vcf.gz
    touch ${prefix}.phased.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        whatshap: \$(whatshap --version)
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
process CUTESV {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cutesv:2.1.1--pyhdfd78af_0' :
        'biocontainers/cutesv:2.1.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(intervals)
    tuple val(meta2), path(fasta), path(fai)
    tuple val(meta3), path(vcf)

    output:
    tuple val(meta), path("${prefix}.vcf"), emit: vcf
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def regions = intervals ? "-include_bed ${intervals}" : ""
    def variants = vcf ? "-Ivcf ${vcf}" : ''

    """
    cuteSV \
        ${bam} \\
        ${fasta} \\
        ${prefix}.vcf \\
        . \\
        --threads $task.cpus \\
        $regions \\
        $variants \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cuteSV: \$( cuteSV --version 2>&1 | sed 's/cuteSV //g' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}" 
    """
    touch "${prefix}.vcf"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cuteSV: \$( cuteSV --version 2>&1 | sed 's/cuteSV //g' )
    END_VERSIONS
    """
}
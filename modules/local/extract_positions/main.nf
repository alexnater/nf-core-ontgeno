process EXTRACT_POSITIONS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.23.0--py39hdd5828d_0' :
        'quay.io/biocontainers/pysam:0.23.0--py39hdd5828d_0' }"

    input:
    tuple val(meta), path(gvcf), path(tbi)
    tuple val(meta2), path(fasta), path(fai)
    tuple val(meta3), path(positions)

    output:
    tuple val(meta), path("*.focal.vcf"), emit: vcf
    path "versions.yml"                 , emit: versions

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    extract_positions_gvcf.py \\
        $gvcf \\
        $fasta \\
        $positions \\
        ${prefix}.focal.gvcf.gz \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
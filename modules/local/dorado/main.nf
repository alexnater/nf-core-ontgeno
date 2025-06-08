process DORADO {
    tag "$meta.id"
    label 'process_gpu'

    container "${projectDir}/assets/containers/dorado_0.9.1.sif"

    input:
    tuple val(meta), path(pod5)
    val(model)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    dorado \\
        basecaller \\
        $args \\
        --emit-fastq \\
        $model \\
        $pod5 | gzip -c > ${prefix}.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(echo \$(dorado --version 2>&1) | tail -n 1)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(echo \$(dorado --version 2>&1) | tail -n 1)
    END_VERSIONS
    """
}
process POD5 {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "containers/pod5_0.3.11.sif"

    input:
    tuple val(meta), path(fast5)

    output:
    tuple val(meta), path("*.pod5"), emit: pod5
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    pod5 \\
        convert fast5 \\
        $args \\
        --threads $task.cpus \\
        --output ${prefix}.pod5 \\
        $fast5

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pod5: \$(echo \$(pod5 --version 2>&1) | sed 's/^.*version: //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pod5

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pod5: \$(echo \$(pod5 --version 2>&1) | sed 's/^.*version: //')
    END_VERSIONS
    """
}
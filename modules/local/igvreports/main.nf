process IGVREPORTS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/igv-reports:1.12.0--pyh7cba7a3_0':
        'biocontainers/igv-reports:1.12.0--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(sites), path(tracks), path(tracks_indices)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.html") , emit: report
    tuple val(meta), path("*.json") , emit: json
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta = fasta ? "--fasta ${fasta}" : ""
    def track_config = tracks.collect { track ->
        if (track.extension == 'bed') {
            "{\n\t\"name\": \"${track.name}\",\n\t\"type\": \"annotation\",\n\t\"format\": \"bed\",\n\t\"url\": \"${track}\"\n}"
        }
        else if (track.extension == 'bam') {
            "{\n\t\"name\": \"${track.name}\",\n\t\"type\": \"alignment\",\n\t\"format\": \"bam\",\n\t\"url\": \"${track}\",\n\t\"indexURL\": \"${track.name}.bai\",\n\t\"showAlignments\": false\n}"
        }
    }.join(',\n')
    def track_arg = tracks ? "--track-config track_config.json" : ""

    """
    echo -e '[\n$track_config\n]' > track_config.json
    
    create_report $sites \\
    $args \\
    $fasta \\
    $track_arg \\
    --output ${prefix}_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        igvreports: \$(python -c "import igv_reports; print(igv_reports.__version__)")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        igvreports: \$(python -c "import igv_reports; print(igv_reports.__version__)")
    END_VERSIONS
    """
}

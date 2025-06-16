process GENOTYPES_TO_VEP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.23.0--py39hdd5828d_0' :
        'quay.io/biocontainers/pysam:0.23.0--py39hdd5828d_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta2), path(tab)
    tuple val(meta3), path(db)
    val(reference)

    output:
    tuple val(meta), path("*.genotypes.tab"), emit: annot
    path "versions.yml"                     , emit: versions

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def db_arg = db ? "--variants $db" : ""
    def ref_arg = reference ? "--reference $reference" : ""

    """
    add_genotypes_to_annotation.py \\
        $tab \\
        $vcf \\
        ${prefix}.genotypes.tab \\
        --samples ${meta.samples.join(' ')} \\
        $db_arg \\
        $ref_arg \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
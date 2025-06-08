process EAGLE2 {
    tag "${meta.id} - ${meta.caller}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${projectDir}/assets/containers/eagle_2.4.1.sif"

    input:
    tuple val(meta) , path(target), path(target_idx)
    tuple val(meta2), path(ref)   , path(ref_idx)
    tuple val(meta3), path(gmap)
    val(chromosome)

    output:
    tuple val(meta), path("${prefix}.bcf")   , emit: bcf, optional: true
    tuple val(meta), path("${prefix}.vcf.gz"), emit: vcf, optional: true
    path "${prefix}.log"                     , emit: log
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.refphased"
    def map_arg = gmap ? "--geneticMapFile ${gmap}" : ''
    if (gmap && args.contains("--geneticMapFile")) {
        error("Genetic map specified both in input channel and ext.args!")
    }
    def chrom_arg = chromosome ? "--chrom ${chromosome}" : ''
    def create_index_target = target.extension == 'bcf' && !target_idx ? 1 : 0
    def create_index_ref = ref.extension == 'bcf' && !ref_idx ? 1 : 0

    """
    if [ $create_index_target -eq 1 ]; then
        bcftools index $target
    fi

    if [ $create_index_ref -eq 1 ]; then
        bcftools index $ref
    fi

    eagle \\
        $args \\
        --numThreads $task.cpus \\
        $map_arg \\
        $chrom_arg \\
        --outPrefix $prefix \\
        --vcfRef $ref \\
        --vcfTarget $target 2>&1 | tee ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eagle2: \$(eagle --help | grep 'Eagle v' | sed -e 's/[| \\t]*Eagle v//' -e 's/[| \\t]\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contain("--vcfOutFormat b") ? "bcf" :
                    args.contain("--vcfOutFormat v") ? "vcf" : "vcf.gz"
    """
    touch ${prefix}.${extension}
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eagle2: \$(eagle --help | grep 'Eagle v' | sed -e 's/[| \\t]*Eagle v//' -e 's/[| \\t]\$//')
    END_VERSIONS
    """
}
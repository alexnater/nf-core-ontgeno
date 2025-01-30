process GATK4_GENOTYPEGVCFS {
    tag "${meta.id} - ${meta.caller}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.6.1.0--py310hdfd78af_0':
        'biocontainers/gatk4:4.6.1.0--py310hdfd78af_0' }"

    input:
    tuple val(meta),  path(input), path(gvcf_index), path(intervals), path(intervals_index)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dict)
    tuple val(meta5), path(dbsnp)
    tuple val(meta6), path(dbsnp_tbi)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.tbi")   , emit: tbi
    path  "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def create_db = input.size() > 1 ? 1 : 0
    def input_command = input.size() == 1 ? "$input" : "gendb://${prefix}"
    def dbsnp_command = dbsnp ? "--dbsnp $dbsnp" : ""
    def interval_command = intervals ? "--intervals $intervals" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK GenotypeGVCFs] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    if [ $create_db -eq 1 ]
    then
        gatk --java-options "-Xmx${avail_mem}M" GenomicsDBImport \\
            ${input.collect(){"--variant $it"}.join(' ')} \\
            --genomicsdb-workspace-path ${prefix} \\
            $interval_command \\
            --tmp-dir . \\
            $args
    fi

    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        GenotypeGVCFs \\
        --variant $input_command \\
        --output ${prefix}.vcf.gz \\
        --reference $fasta \\
        $interval_command \\
        $dbsnp_command \\
        --tmp-dir . \\
        $args2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}

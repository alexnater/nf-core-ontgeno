/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { NANOPLOT               } from '../modules/nf-core/nanoplot'
include { CAT_CAT                } from '../modules/nf-core/cat/cat'
include { PREPROCESSING          } from '../subworkflows/local/preprocessing'
include { QC                     } from '../subworkflows/local/qc'
include { QC as QC_FILTERED      } from '../subworkflows/local/qc'
include { MINIMAP2_ALIGN         } from '../modules/nf-core/minimap2/align'
include { SAMTOOLS_VIEW          } from '../modules/nf-core/samtools/view'
include { MULTIQC                } from '../modules/nf-core/multiqc'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { BASECALLING            } from '../subworkflows/local/basecalling'
include { BAM_STATS              } from '../subworkflows/local/bam_stats'
include { VARIANT_CALLING        } from '../subworkflows/local/variant_calling'
include { PHASE_VARIANTS         } from '../subworkflows/local/phase_variants_single'
include { ANNOTATE_SNPS          } from '../subworkflows/local/annotate_snps'
include { SV_CALLING             } from '../subworkflows/local/svcalling'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_ontgeno_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ONTGENO {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    // define variables
    def ref_file = file(params.fasta, checkIfExists: true)
    def fai_file = file("${params.fasta}.fai", checkIfExists: true)
    def dict_file = file(ref_file.parent / "${ref_file.baseName}.dict", checkIfExists: true)
    def bed_file = file(params.bed, checkIfExists: true)
    def str_file = params.str_file ? file(params.str_file, checkIfExists: true) : []
    def model_file = file(params.genotype_model, checkIfExists: true)
    def config_file = file(params.glnexus_config, checkIfExists: true)
    def pos_file = params.focal_positions ? file(params.focal_positions, checkIfExists: true) : null
    def panel_file = params.panel ? file(params.panel, checkIfExists: true) : null
    def vep_cache = params.vep_cache ? file(params.vep_cache, type: 'dir', checkIfExists: true) : null
    def ch_fasta_fai = channel.value([ [id:'ref'], ref_file, fai_file ])
    def ch_dict = channel.value([ [id:'ref'], dict_file ])

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    ch_samplesheet.branch { meta, infiles ->
        fast5:  meta.fast5
        bam:    meta.bam
        fq_dir: meta.fq_dir
        fastq:  true
    }.set { ch_input }

    //
    // MODULE: Run fastq concatenation
    //
    CAT_CAT (
        ch_input.fq_dir
    )
    ch_versions = ch_versions.mix(CAT_CAT.out.versions.first())

    //
    // SUBWORKFLOW: basecalling
    //
    BASECALLING (
        ch_input.fast5,
        params.basecalling_model
    ).fastq
     .mix(ch_input.fastq, CAT_CAT.out.file_out)
     .set { ch_fastq }
    ch_versions = ch_versions.mix(BASECALLING.out.versions)

    //
    // SUBWORKFLOW: qc
    //
    QC (
        ch_fastq
    )
    ch_multiqc_files = ch_multiqc_files.mix(QC.out.fastqc.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(QC.out.nanoplot.collect{it[1]})
    ch_versions = ch_versions.mix(QC.out.versions)

    //
    // SUBWORKFLOW: preprocessing
    //
    PREPROCESSING (
        ch_fastq,
        params.trimmer
    )
    ch_multiqc_files = ch_multiqc_files.mix(PREPROCESSING.out.json.collect{it[1]})
    ch_versions = ch_versions.mix(PREPROCESSING.out.versions)

    //
    // SUBWORKFLOW: qc
    //
    QC_FILTERED (
        PREPROCESSING.out.trimmed
    )
    ch_versions = ch_versions.mix(QC_FILTERED.out.versions)

    //
    // MODULE: Run minimap2
    //
    PREPROCESSING.out.trimmed
        .map { meta, fastq -> [ [id: meta.sample, sample: meta.sample, model: meta.model], fastq ] }
        .groupTuple()
        .set { ch_fastqs }
    
    MINIMAP2_ALIGN (
        ch_fastqs,
        ch_fasta_fai.map { meta, fasta, fai -> [ meta, fasta ] },
        true,
        "bai",
        false,
        false        
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())

    MINIMAP2_ALIGN.out.bam
        .join(MINIMAP2_ALIGN.out.index, failOnDuplicate:true, failOnMismatch:true)
        .set { ch_mapped }

    //
    // SUBWORKFLOW: bam_stats
    //
    BAM_STATS (
        ch_mapped,
        ch_fasta_fai,
        PREPROCESSING.out.json,
        bed_file
    )
    ch_versions = ch_versions.mix(BAM_STATS.out.versions)

    //
    // MODULE: Run samtools view
    //
    SAMTOOLS_VIEW (
        ch_mapped,
        ch_fasta_fai.map { meta, fasta, fai -> [ meta, fasta ] },
        [],
        "bai"
    )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions.first())

    SAMTOOLS_VIEW.out.bam
        .join(SAMTOOLS_VIEW.out.bai, failOnDuplicate:true, failOnMismatch:true)
        .set { ch_bam_bai }

    //
    // SUBWORKFLOW: variant_calling
    //
    VARIANT_CALLING (
        ch_bam_bai,
        ch_fasta_fai,
        ch_dict,
        bed_file,
        params.genotype_caller,
        model_file,
        config_file,
        pos_file
    )
    ch_multiqc_files = ch_multiqc_files.mix(VARIANT_CALLING.out.stats.collect{it[1]})
    ch_versions = ch_versions.mix(VARIANT_CALLING.out.versions)

    //
    // SUBWORKFLOW: phase_variants
    //
    PHASE_VARIANTS (
        VARIANT_CALLING.out.vcf_tbi_single,
        ch_bam_bai,
        ch_fasta_fai,
        bed_file,
        panel_file
    )
    ch_multiqc_files = ch_multiqc_files.mix(PHASE_VARIANTS.out.stats.collect{it[1]})
    ch_versions = ch_versions.mix(PHASE_VARIANTS.out.versions)

    //
    // SUBWORKFLOW: annotate_snps
    //
    ANNOTATE_SNPS (
        PHASE_VARIANTS.out.vcf_tbi,
        ch_fasta_fai,
        params.genome,
        params.species,
        params.vep_cache_version,
        vep_cache
    )
    ch_multiqc_files = ch_multiqc_files.mix(ANNOTATE_SNPS.out.report.collect{it[1]})
    ch_versions = ch_versions.mix(ANNOTATE_SNPS.out.versions)

    //
    // SUBWORKFLOW: svcalling
    //
    SV_CALLING (
        ch_bam_bai,
        ch_fasta_fai,
        bed_file,
        str_file
    )
    ch_versions = ch_versions.mix(SV_CALLING.out.versions)

    //
    // Collate and save software versions
    //
    def topic_versions = Channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'ontgeno_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

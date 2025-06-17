/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { NANOPLOT               } from '../modules/nf-core/nanoplot/main'
include { FASTP                  } from '../modules/nf-core/fastp/main'
include { MINIMAP2_ALIGN         } from '../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_VIEW          } from '../modules/nf-core/samtools/view'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { BASECALLING            } from '../subworkflows/local/basecalling'
include { BAM_STATS              } from '../subworkflows/local/bam_stats'
include { VARIANT_CALLING        } from '../subworkflows/local/variant_calling'
include { PHASE_VARIANTS         } from '../subworkflows/local/phase_variants'
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

def ref_file = file(params.fasta, checkIfExists: true)
def fai_file = file("${params.fasta}.fai", checkIfExists: true)
def dict_file = file(ref_file.parent / "${ref_file.baseName}.dict", checkIfExists: true)
def bed_file = file(params.bed, checkIfExists: true)
def str_file = params.str_file ? file(params.str_file, checkIfExists: true) : []
def model_file = file(params.genotype_model, checkIfExists: true)
def config_file = file(params.glnexus_config, checkIfExists: true)
def pos_file = file(params.focal_positions, checkIfExists: true)
def panel_file = params.panel ? file(params.panel, checkIfExists: true) : null
def vep_cache = params.vep_cache ? file(params.vep_cache, type: 'dir', checkIfExists: true) : null
def ch_fasta_fai = Channel.value([ [id:'ref'], ref_file, fai_file ])
def ch_dict = Channel.value([ [id:'ref'], dict_file ])


workflow ONTGENO {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_samplesheet.branch { meta, infiles ->
        fast5:  meta.fast5
        bam:    meta.bam
        fq_dir: meta.fq_dir
        fastq:  true
    }.set { ch_input }

    ch_input.fq_dir
        .map { meta, fastq_files ->
            [ meta, 1..fastq_files.size(), fastq_files ]
        }
        .transpose()
        .map { meta, rep, fastq -> [ meta + [id: "${meta.id}_${rep}"], fastq ] }
        .set { ch_fastq_dir }

    //
    // SUBWORKFLOW: basecalling
    //
    BASECALLING (
        ch_input.fast5,
        params.basecalling_model
    ).fastq
     .mix(ch_input.fastq, ch_fastq_dir)
     .set { ch_fastq }
    ch_versions = ch_versions.mix(BASECALLING.out.versions)

    //
    // MODULE: Run fastp
    //
    FASTP (
        ch_fastq,
        [],
        false,
        false,
        false
    )
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        FASTP.out.reads
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: Run NanoPlot
    //
    NANOPLOT (
        FASTP.out.reads
    )
    ch_multiqc_files = ch_multiqc_files.mix(NANOPLOT.out.txt.collect{it[1]})
    ch_versions = ch_versions.mix(NANOPLOT.out.versions.first())

    //
    // MODULE: Run minimap2
    //
    FASTP.out.reads
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
        FASTP.out.json,
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
        model_file,
        config_file,
        pos_file
    )
    ch_versions = ch_versions.mix(VARIANT_CALLING.out.versions)

    //
    // SUBWORKFLOW: phase_variants
    //
    PHASE_VARIANTS (
        VARIANT_CALLING.out.vcf_tbi,
        VARIANT_CALLING.out.bam_bai,
        ch_fasta_fai,
        bed_file,
        panel_file
    )
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
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'ontgeno_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
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

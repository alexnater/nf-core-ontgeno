/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EAGLE2                 } from '../../../modules/local/eagle2'
include { WHATSHAP_PHASE         } from '../../../modules/local/whatshap/phase'
include { WHATSHAP_STATS         } from '../../../modules/local/whatshap/stats'
include { IGVREPORTS             } from '../../../modules/nf-core/igvreports'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PHASE_VARIANTS {

    take:
    ch_vcf_tbi    // channel: [ meta, vcf, tbi ]
    ch_bam_bai    // channel: [ meta, [bam], [bai] ]
    ch_fasta_fai  // channel: [ meta, fasta, fai ]
    bed_file      // BED file with genomic intervals
    panel_file    // haplotype reference panel file

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Run Eagle2
    //
    EAGLE2 (
        ch_vcf_tbi,
        [ [id: 'ref'], panel_file, [] ],
        [ [:], [] ],
        ""
    )
    ch_versions = ch_versions.mix(EAGLE2.out.versions.first())

    //
    // MODULE: Run WhatsHap phase
    //
    ch_vcf_tbi
        .join(EAGLE2.out.vcf, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, vcf, tbi, phased_vcf -> [ meta - meta.subMap('typer'), meta, vcf, tbi, phased_vcf ] }
        .combine(ch_bam_bai, by: 0)
        .multiMap { groupkey, meta, vcf, tbi, phased_vcf, bams, bais ->
            calls: [ meta, vcf ]
            reads: [ meta, [*bams, phased_vcf], bais ]
        }.set { ch_to_phase }
  
    WHATSHAP_PHASE (
        ch_to_phase.calls,
        ch_to_phase.reads,
        ch_fasta_fai
    )
    ch_versions = ch_versions.mix(WHATSHAP_PHASE.out.versions.first())

    //
    // MODULE: Run WhatsHap stats
    //
    WHATSHAP_STATS (
        WHATSHAP_PHASE.out.vcf_tbi
            .map { meta, vcf, tbi -> [ meta, vcf, tbi, meta.samples ] }
    )
    ch_versions = ch_versions.mix(WHATSHAP_STATS.out.versions.first())

    //
    // MODULE: Run igv-reports
    //
    WHATSHAP_PHASE.out.vcf_tbi
        .map { meta, vcf, tbi -> [ meta.subMap(['id', 'caller', 'typer']), vcf, tbi ] }
        .join(WHATSHAP_STATS.out.bed
                .map { meta, bed -> [ meta.subMap(['id', 'caller', 'typer']), bed ] },
              failOnDuplicate:true, failOnMismatch:true)
        .map { meta, vcf, tbi, bed -> [ meta, vcf, [bed], [tbi] ] }
        .set { ch_to_reports }

    IGVREPORTS (
        ch_to_reports,
        ch_fasta_fai
    )
    ch_versions = ch_versions.mix(IGVREPORTS.out.versions.first())


    emit:
    vcf_tbi  = WHATSHAP_PHASE.out.vcf_tbi   // channel: [ meta, vcf, tbi ]
    gtf      = WHATSHAP_STATS.out.gtf       // channel: [ meta, gtf ]
    report   = IGVREPORTS.out.report        // channel: [ meta, html ]
    versions = ch_versions                  // channel: [ path(versions.yml) ]
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SNIFFLES as SNIFFLES_SINGLE       } from '../../../modules/nf-core/sniffles/main'
include { SNIFFLES as SNIFFLES_JOINT        } from '../../../modules/nf-core/sniffles/main'
include { CUTESV as CUTESV_SINGLE           } from '../../../modules/nf-core/cutesv/main'
include { CUTESV as CUTESV_GENOTYPE         } from '../../../modules/nf-core/cutesv/main'
include { SURVIVOR_MERGE                    } from '../../../modules/nf-core/survivor/merge/main'
include { SURVIVOR_MERGE as SURVIVOR_MERGE2 } from '../../../modules/nf-core/survivor/merge/main'
include { IGVREPORTS as IGVREPORTS_SNIFFLES } from '../../../modules/nf-core/igvreports/main'
include { IGVREPORTS as IGVREPORTS_CUTESV   } from '../../../modules/nf-core/igvreports/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SV_CALLING {

    take:
    ch_bam_bai         // channel: [ meta, bam, bai ]
    ch_fasta_fai       // channel: [ meta, fasta, fai ]
    bed_file           // BED file with genomic intervals
    tandem_repeats     // BED file with tandem repeat annotations

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Run Sniffles per individual
    //
    SNIFFLES_SINGLE (
        ch_bam_bai,
        ch_fasta_fai.map { meta, fasta, fai -> [ meta, fasta ] },
        [ [id:'repeats'], tandem_repeats ],
        false,
        true
    ).snf
     .map { meta, snf -> [ [id: 'joint', caller: 'sniffles'], snf ] }
     .groupTuple(sort: true)
     .set { ch_snf }
    ch_versions = ch_versions.mix(SNIFFLES_SINGLE.out.versions.first())

    //
    // MODULE: Run Sniffles jointly
    //
    SNIFFLES_JOINT (
        ch_snf.map { meta, snf -> [ meta, snf, [] ] },
        ch_fasta_fai.map { meta, fasta, fai -> [ meta, fasta ] },
        [ [id:'repeats'], tandem_repeats ],
        true,
        false
    )
    ch_versions = ch_versions.mix(SNIFFLES_JOINT.out.versions.first())

    //
    // MODULE: Run igv-reports
    //
    IGVREPORTS_SNIFFLES (
        SNIFFLES_JOINT.out.vcf
            .map { meta, vcf -> [ meta, vcf, [], [] ] },
        ch_fasta_fai
    )
    ch_versions = ch_versions.mix(IGVREPORTS_SNIFFLES.out.versions.first())

    //
    // MODULE: Run cuteSV per individual in discovery mode
    //
    CUTESV_SINGLE (
        ch_bam_bai.map { meta, bam, bai ->
            [ meta, bam, bai, bed_file ]
        },
        ch_fasta_fai,
        [ [:], [] ]
    ).vcf
     .map { meta, vcf -> [ [id: 'joint', caller: 'cutesv'], vcf ] }
     .groupTuple(sort: true)
     .set { ch_vcf }
    ch_versions = ch_versions.mix(CUTESV_SINGLE.out.versions.first())

    //
    // MODULE: Merge VCFs with SURVIVOR
    //
    SURVIVOR_MERGE (
        ch_vcf,
        1000,
        1,
        1,
        1,
        1,
        2
    )
    ch_versions = ch_versions.mix(SURVIVOR_MERGE.out.versions.first())

    //
    // MODULE: Run cuteSV per individual in genotype mode
    //
    CUTESV_GENOTYPE (
        ch_bam_bai.map { meta, bam, bai ->
            [ meta, bam, bai, bed_file ]
        },
        ch_fasta_fai,
        SURVIVOR_MERGE.out.vcf.first()
    ).vcf
     .map { meta, vcf -> [ [id: 'genotyped', caller: 'cutesv'], vcf ] }
     .groupTuple(sort: true)
     .set { ch_genotyped }
    ch_versions = ch_versions.mix(CUTESV_GENOTYPE.out.versions.first())

    //
    // MODULE: Merge genotyped VCFs with SURVIVOR
    //
    SURVIVOR_MERGE2 (
        ch_genotyped,
        1000,
        1,
        1,
        1,
        1,
        2
    )
    ch_versions = ch_versions.mix(SURVIVOR_MERGE2.out.versions.first())

    //
    // MODULE: Run igv-reports
    //
    IGVREPORTS_CUTESV (
        SURVIVOR_MERGE2.out.vcf
            .map { meta, vcf -> [ meta, vcf, [], [] ] },
        ch_fasta_fai
    )
    ch_versions = ch_versions.mix(IGVREPORTS_CUTESV.out.versions.first())


    emit:
    vcf      = SNIFFLES_JOINT.out.vcf           // channel: [ meta, vcf ]
    report   = IGVREPORTS_SNIFFLES.out.report   // channel: [ meta, html ]
    versions = ch_versions                      // channel: [ path(versions.yml) ]
}
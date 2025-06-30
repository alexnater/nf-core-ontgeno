/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CLAIR3                               } from '../../../modules/local/clair3'
include { TABIX_TABIX                          } from '../../../modules/nf-core/tabix/tabix'
include { EXTRACT_POSITIONS as FOCALPOS_CLAIR3 } from '../../../modules/local/extract_positions'
include { DEEPVARIANT_RUNDEEPVARIANT           } from '../../../modules/nf-core/deepvariant/rundeepvariant'
include { EXTRACT_POSITIONS as FOCALPOS_DV     } from '../../../modules/local/extract_positions'
include { GATK4_HAPLOTYPECALLER                } from '../../../modules/nf-core/gatk4/haplotypecaller'
include { GATK4_COMBINEGVCFS                   } from '../../../modules/nf-core/gatk4/combinegvcfs'
include { GATK4_GENOTYPEGVCFS                  } from '../../../modules/nf-core/gatk4/genotypegvcfs'
include { GLNEXUS                              } from '../../../modules/nf-core/glnexus'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VARIANT_CALLING {

    take:
    ch_bam_bai    // channel: [ meta, bam, bai ]
    ch_fasta_fai  // channel: [ meta, fasta, fai ]
    ch_dict       // channel: [ meta, dict ]
    bed_file      // BED file with genomic intervals
    model_file    // Clair3 model file
    config_file   // GL Nexus config file
    pos_file      // File with focal positions

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Run Clair3
    //
    CLAIR3 (
        ch_bam_bai.map { meta, bam, bai ->
            [ meta, bam, bai, bed_file ]
        },
        ch_fasta_fai,
        "ont",
        ch_bam_bai.map { meta, bam, bai -> meta.model ?: model_file }
    )
    ch_versions = ch_versions.mix(CLAIR3.out.versions.first())

    //
    // MODULE: Run Tabix
    //
    TABIX_TABIX (
        CLAIR3.out.gvcf
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    // Join with indices
    CLAIR3.out.gvcf
        .join(TABIX_TABIX.out.tbi, failOnDuplicate:true, failOnMismatch:true)
        .set { ch_gvcf_tbi_clair3 }

    //
    // MODULE: Run extract_positions
    //
    if (pos_file) {
        FOCALPOS_CLAIR3 (
            ch_gvcf_tbi_clair3,
            ch_fasta_fai,
            [ [:], pos_file ]
        )
        ch_versions = ch_versions.mix(FOCALPOS_CLAIR3.out.versions.first())
    }

    // Join with bam files and group
    ch_gvcf_tbi_clair3
        .join(ch_bam_bai, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, gvcf, tbi, bam, bai -> [ [id: 'joint', caller: 'clair3'], [ meta.sample, gvcf, tbi, bam, bai ] ] }
        .groupTuple(sort: { a, b -> a[0] <=> b[0] })
        .map { meta, tuples ->
            def (samples, gvcfs, tbis, bams, bais) = tuples.transpose()
            [ meta + [samples: tuple(samples)], gvcfs, tbis, bams, bais ]
        }
        .set { ch_from_clair3 }

    //
    // MODULE: Run DeepVariant
    //
    DEEPVARIANT_RUNDEEPVARIANT (
        ch_bam_bai.map { meta, bam, bai ->
            [ meta, bam, bai, bed_file ]
        },
        ch_fasta_fai.map { meta, fasta, fai -> [ meta, fasta ] },
        ch_fasta_fai.map { meta, fasta, fai -> [ meta, fai ] },
        [ [id:'ref'], [] ],
        [ [:], [] ]
    )
    ch_versions = ch_versions.mix(DEEPVARIANT_RUNDEEPVARIANT.out.versions.first())

    // Join with indices
    DEEPVARIANT_RUNDEEPVARIANT.out.gvcf
        .join(DEEPVARIANT_RUNDEEPVARIANT.out.gvcf_tbi, failOnDuplicate:true, failOnMismatch:true)
        .set { ch_gvcf_tbi_dv }

    //
    // MODULE: Run extract_positions
    //
    if (pos_file) {
        FOCALPOS_DV (
            ch_gvcf_tbi_dv,
            ch_fasta_fai,
            [ [:], pos_file ]
        )
        ch_versions = ch_versions.mix(FOCALPOS_DV.out.versions.first())
    }

    // Join with bam files and group
    ch_gvcf_tbi_dv
        .join(ch_bam_bai, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, gvcf, tbi, bam, bai -> [ [id: 'joint', caller: 'deepvariant'], [ meta.sample, gvcf, tbi, bam, bai ] ] }
        .groupTuple(sort: { a, b -> a[0] <=> b[0] })
        .map { meta, tuples ->
            def (samples, gvcfs, tbis, bams, bais) = tuples.transpose()
            [ meta + [samples: tuple(samples)], gvcfs, tbis, bams, bais ]
        }
        .set { ch_from_dv }

    //
    // MODULE: Run GATK HaplotypeCaller
    //
    Channel.empty().set { ch_from_gatk }

/*
    GATK4_HAPLOTYPECALLER (
        ch_bam_bai.map { meta, bam, bai -> [ meta, bam, bai, bed_file, [] ] },
        ch_fasta_fai.map { meta, fasta, fai -> [ meta, fasta ] },
        ch_fasta_fai.map { meta, fasta, fai -> [ meta, fai ] },
        ch_dict,
        [ [:], [] ],
        [ [:], [] ],
    )
    ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions.first())

    // Join with indices and group
    GATK4_HAPLOTYPECALLER.out.vcf
        .join(GATK4_HAPLOTYPECALLER.out.tbi, failOnDuplicate:true, failOnMismatch:true)
        .join(ch_bam_bai, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, gvcf, tbi, bam, bai -> [ [id: 'joint', caller: 'haplotypecaller'], [ meta.sample, gvcf, tbi, bam, bai ] ] }
        .groupTuple(sort: { a, b -> a[0] <=> b[0] })
        .map { meta, tuples ->
            def (samples, gvcfs, tbis, bams, bais) = tuples.transpose()
            [ meta + [samples: tuple(samples)], gvcfs, tbis, bams, bais ]
        }
        .set { ch_from_gatk }
*/

    ch_from_clair3
        .mix(ch_from_dv, ch_from_gatk)
        .multiMap { meta, gvcfs, tbis, bams, bais ->
            to_gatk:    [ meta, gvcfs, tbis ]
            to_glnexus: [ meta, gvcfs ]
            bam_bai:    [ meta, bams, bais ]
        }.set { ch_calls }

    //
    // MODULE: Run GATK CombineGVCFs
    //
    Channel.empty().set { ch_vcf_tbi }

/*
    GATK4_COMBINEGVCFS (
        ch_calls.to_gatk.filter { meta, gvcfs, tbis -> meta.caller != 'deepvariant' },
        ch_fasta_fai.map { meta, fasta, fai -> fasta },
        ch_fasta_fai.map { meta, fasta, fai -> fai },
        ch_dict.map { meta, dict -> dict }
    )
    ch_versions = ch_versions.mix(GATK4_COMBINEGVCFS.out.versions.first())

    GATK4_COMBINEGVCFS.out.combined_gvcf
        .join(GATK4_COMBINEGVCFS.out.combined_tbi, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, gvcf, tbi -> [ meta, gvcf, tbi, bed_file, [] ] }
        .set { ch_to_genotype }

    //
    // MODULE: Run GATK GenotypeGVCFs
    //
    GATK4_GENOTYPEGVCFS (
        ch_to_genotype,
        ch_fasta_fai.map { meta, fasta, fai -> [ meta, fasta ] },
        ch_fasta_fai.map { meta, fasta, fai -> [ meta, fai ] },
        ch_dict,
        [ [:], [] ],
        [ [:], [] ]
    )
    ch_versions = ch_versions.mix(GATK4_GENOTYPEGVCFS.out.versions.first())

    GATK4_GENOTYPEGVCFS.out.vcf
        .join(GATK4_GENOTYPEGVCFS.out.tbi, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, vcf, tbi -> [ meta + [typer: 'gatk'], vcf, tbi ] }
        .set { ch_vcf_tbi }
*/
    
    // Make sure to use the correct config for GLnexus
    ch_calls.to_glnexus
        .multiMap { meta, gvcfs ->
            input:  [ meta, gvcfs ]
            preset: meta.caller == 'clair3' ? null : (meta.caller == 'deepvariant' ? 'DeepVariant' : 'gatk')
            config: meta.caller == 'clair3' ? config_file : []
        }.set { ch_glnexus }

    //
    // MODULE: Run GLnexus
    //
    GLNEXUS (
        ch_glnexus.input,
        [ [:], bed_file ],
        ch_glnexus.preset,
        ch_glnexus.config
    )
    ch_versions = ch_versions.mix(GLNEXUS.out.versions.first())

    ch_vcf_tbi
        .mix(GLNEXUS.out.bcf.map { meta, bcf -> [ meta + [typer: 'glnexus'], bcf, [] ] })
        .set { vcf_tbi }


    emit:
    vcf_tbi                                 // channel: [ meta, vcf, tbi ]
    bam_bai  = ch_calls.bam_bai             // channel: [ meta, [bam], [bai] ]
    versions = ch_versions                  // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CLAIR3                               } from '../../../modules/local/clair3'
include { TABIX_TABIX as TABIX_VCF             } from '../../../modules/nf-core/tabix/tabix'
include { TABIX_TABIX as TABIX_GVCF            } from '../../../modules/nf-core/tabix/tabix'
include { DEEPVARIANT_RUNDEEPVARIANT           } from '../../../modules/nf-core/deepvariant/rundeepvariant'
include { EXTRACT_POSITIONS                    } from '../../../modules/local/extract_positions'
include { BCFTOOLS_MPILEUP as BCFTOOLS_JOINT   } from '../../../modules/local/bcftools'
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
    caller        // Genotype caller to use
    model_file    // Clair3 model file
    config_file   // GL Nexus config file
    pos_file      // File with focal positions

    main:

    ch_versions = channel.empty()
    ch_gvcf_tbi = channel.empty()
    ch_stats = channel.empty()

    if (caller == 'clair3') {
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
        // MODULE: Run Tabix on gvcf output
        //
        TABIX_GVCF (
            CLAIR3.out.gvcf
        )
        ch_versions = ch_versions.mix(TABIX_GVCF.out.versions.first())

        // Join with indices
        CLAIR3.out.gvcf
            .join(TABIX_GVCF.out.tbi, failOnDuplicate:true, failOnMismatch:true)
            .set { ch_gvcf_tbi }

    } else if (caller == 'deepvariant') {
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

        // Join gvcf output with indices
        DEEPVARIANT_RUNDEEPVARIANT.out.gvcf
            .join(DEEPVARIANT_RUNDEEPVARIANT.out.gvcf_tbi, failOnDuplicate:true, failOnMismatch:true)
            .set { ch_gvcf_tbi }
    
    } else if (caller != 'bcftools') {
        error "Invalid genotype caller specified: ${caller}"
    }

    //
    // MODULE: Run extract_positions
    //
    if (pos_file) {
        EXTRACT_POSITIONS (
            ch_gvcf_tbi,
            ch_fasta_fai,
            [ [:], pos_file ]
        )
        ch_versions = ch_versions.mix(EXTRACT_POSITIONS.out.versions.first())
    }

    // Run joint variant calling
    if (caller == 'clair3' || caller == 'deepvariant') {

        // Group and sort gvcf files and make sure to use the correct config for GLnexus
        def ch_glnexus = ch_gvcf_tbi
            .map { meta, gvcf, tbi -> [ [id: 'joint'], [ meta.sample, gvcf, tbi ] ] }
            .groupTuple(sort: { a, b -> a[0] <=> b[0] })
            .multiMap { meta, tuples ->
                def (samples, gvcfs, tbis) = tuples.transpose()
                input: [ meta + [samples: tuple(samples)], gvcfs ]
                preset: caller == 'clair3' ? null : 'DeepVariant'
                config: caller == 'clair3' ? config_file : []
            }
        
        //
        // MODULE: Run GLnexus
        //
        GLNEXUS (
            ch_glnexus.input,
            [ [:], bed_file ],
            ch_glnexus.preset,
            ch_glnexus.config
        ).bcf
            .map { meta, bcf -> [ meta, bcf, [] ] }
            .set { vcf_tbi }
        ch_versions = ch_versions.mix(GLNEXUS.out.versions.first())

    } else {
        // Group and sort bam files
        def ch_to_bcftools = ch_bam_bai
            .map { meta, bam, bai -> [ [id: 'joint'], [ meta.sample, bam, bai ] ] }
            .groupTuple(sort: { a, b -> a[0] <=> b[0] })
            .map { meta, tuples ->
                def (samples, bams, bais) = tuples.transpose()
                [ meta + [samples: tuple(samples)], bams, bais, bed_file ]
            }

        //
        // MODULE: Run bcftools mpileup over all samples
        //
        BCFTOOLS_JOINT (
            ch_to_bcftools,
            ch_fasta_fai.map { meta, fasta, fai -> [ meta, fasta ] }
        )
        .vcf_tbi
        .set { vcf_tbi }

        ch_stats.mix(BCFTOOLS_JOINT.out.stats)
    }

    emit:
    vcf_tbi                           // channel: [ meta, vcf, tbi ]
    gvcf_tbi = ch_gvcf_tbi            // channel: [ meta, gvcf, tbi ]
    stats = ch_stats                  // channel: [ meta, txt ]
    versions = ch_versions            // channel: [ path(versions.yml) ]
}

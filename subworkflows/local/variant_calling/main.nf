/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CLAIR3                               } from '../../../modules/local/clair3'
include { TABIX_TABIX as TABIX_VCF             } from '../../../modules/nf-core/tabix/tabix'
include { TABIX_TABIX as TABIX_GVCF            } from '../../../modules/nf-core/tabix/tabix'
include { EXTRACT_POSITIONS as FOCALPOS_CLAIR3 } from '../../../modules/local/extract_positions'
include { DEEPVARIANT_RUNDEEPVARIANT           } from '../../../modules/nf-core/deepvariant/rundeepvariant'
include { EXTRACT_POSITIONS as FOCALPOS_DV     } from '../../../modules/local/extract_positions'
include { GATK4_HAPLOTYPECALLER                } from '../../../modules/nf-core/gatk4/haplotypecaller'
include { GATK4_COMBINEGVCFS                   } from '../../../modules/nf-core/gatk4/combinegvcfs'
include { GATK4_GENOTYPEGVCFS                  } from '../../../modules/nf-core/gatk4/genotypegvcfs'
include { BCFTOOLS_MPILEUP                     } from '../../../modules/nf-core/bcftools/mpileup'
include { BCFTOOLS_MPILEUP as BCFTOOLS_JOINT   } from '../../../modules/nf-core/bcftools/mpileup'
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
        // MODULE: Run Tabix on vcf output
        //
        TABIX_VCF (
            CLAIR3.out.vcf
        )
        ch_versions = ch_versions.mix(TABIX_VCF.out.versions.first())

        // Join with indices
        CLAIR3.out.vcf
            .join(TABIX_VCF.out.tbi, failOnDuplicate:true, failOnMismatch:true)
            .set { vcf_tbi_single }

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

        //
        // MODULE: Run extract_positions
        //
        if (pos_file) {
            FOCALPOS_CLAIR3 (
                ch_gvcf_tbi,
                ch_fasta_fai,
                [ [:], pos_file ]
            )
            ch_versions = ch_versions.mix(FOCALPOS_CLAIR3.out.versions.first())
        }

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

        // Join vcf output with indices
        DEEPVARIANT_RUNDEEPVARIANT.out.vcf
            .join(DEEPVARIANT_RUNDEEPVARIANT.out.vcf_tbi, failOnDuplicate:true, failOnMismatch:true)
            .set { vcf_tbi_single }

        // Join gvcf output with indices
        DEEPVARIANT_RUNDEEPVARIANT.out.gvcf
            .join(DEEPVARIANT_RUNDEEPVARIANT.out.gvcf_tbi, failOnDuplicate:true, failOnMismatch:true)
            .set { ch_gvcf_tbi }

        //
        // MODULE: Run extract_positions
        //
        if (pos_file) {
            FOCALPOS_DV (
                ch_gvcf_tbi,
                ch_fasta_fai,
                [ [:], pos_file ]
            )
            ch_versions = ch_versions.mix(FOCALPOS_DV.out.versions.first())
        }
    
    } else if (caller == 'bcftools') {
        //
        // MODULE: Run bcftools mpileup per sample
        //
        BCFTOOLS_MPILEUP (
            ch_bam_bai.map { meta, bam, bai -> [ meta, bam, bed_file ] },
            ch_fasta_fai.map { meta, fasta, fai -> [ meta, fasta ] },
            false
        )

        // Join vcf with indices
        BCFTOOLS_MPILEUP.out.vcf
            .join(BCFTOOLS_MPILEUP.out.tbi, failOnDuplicate:true, failOnMismatch:true)
            .set { vcf_tbi_single }

        ch_stats.mix(BCFTOOLS_MPILEUP.out.stats)

    } else {
        error "Invalid genotype caller specified: ${caller}"
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
            .set { vcf_tbi_joint }
        ch_versions = ch_versions.mix(GLNEXUS.out.versions.first())

    } else {
        channel.empty()
            .set{ vcf_tbi_joint } 
/*
        // Group and sort bam files:
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
            ch_fasta_fai.map { meta, fasta, fai -> [ meta, fasta ] },
            false
        )

        // Join vcf with indices
        BCFTOOLS_JOINT.out.vcf
            .join(BCFTOOLS_JOINT.out.tbi, failOnDuplicate:true, failOnMismatch:true)
            .set { vcf_tbi_joint }

        ch_stats.mix(BCFTOOLS_JOINT.out.stats)
*/
    }

    emit:
    vcf_tbi_single                          // channel: [ meta, vcf, tbi ]
    vcf_tbi_joint                           // channel: [ meta, vcf, tbi ]
    stats = ch_stats                        // channel: [ meta, txt ]
    versions = ch_versions                  // channel: [ path(versions.yml) ]
}

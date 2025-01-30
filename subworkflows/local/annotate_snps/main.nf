/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ENSEMBLVEP_DOWNLOAD } from '../../../modules/nf-core/ensemblvep/download/main'
include { ENSEMBLVEP_VEP      } from '../../../modules/nf-core/ensemblvep/vep/main'
include { TABIX_TABIX         } from '../../../modules/nf-core/tabix/tabix/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ANNOTATE_SNPS {

    take:
    ch_vcf_tbi     // channel: [ meta, vcf, tbi ]
    ch_fasta_fai   // value channel: [ meta, fasta, fai ]
    genome         // Genome version
    species        // Species to annotate
    cache_version  // Ensembl version of cache
    cache          // optional: path to vep cache

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: ENSEMBLVEP_DOWNLOAD if needed
    //
    if (cache) {
        Channel.value(cache).set { ch_cache }
    } else {
        ENSEMBLVEP_DOWNLOAD (
            [ [id:'cache'], genome, species, cache_version ]
        ).cache
         .map { meta, cache -> cache }
         .first()
         .set { ch_cache }
        ch_versions = ch_versions.mix(ENSEMBLVEP_DOWNLOAD.out.versions.first())
    }

    //
    // MODULE: ENSEMBLVEP_VEP
    //
    ENSEMBLVEP_VEP (
        ch_vcf_tbi.map { meta, vcf, tbi -> [ meta, vcf, [] ] },
        genome,
        species,
        cache_version,
        ch_cache,
        ch_fasta_fai.map { meta, fasta, fai -> [ meta, fasta ] },
        []
    )
    ch_versions = ch_versions.mix(ENSEMBLVEP_VEP.out.versions.first())


    emit:
    vcf      = ENSEMBLVEP_VEP.out.vcf       // channel: [ meta, vcf ]
    tab      = ENSEMBLVEP_VEP.out.tab       // channel: [ meta, tab ]
    json     = ENSEMBLVEP_VEP.out.json      // channel: [ meta, json ]
    report   = ENSEMBLVEP_VEP.out.report    // channel: html
    versions = ch_versions                  // channel: [ path(versions.yml) ]
}
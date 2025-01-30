/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { POD5   } from '../../../modules/local/pod5/main'
include { DORADO } from '../../../modules/local/dorado/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BASECALLING {

    take:
    ch_fast5       // channel: [ meta, fast5s ]
    model          // Dorado model

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Run POD5 conversion
    //
    POD5(ch_fast5)
    ch_versions = ch_versions.mix(POD5.out.versions.first())

    //
    // MODULE: Run Dorado
    //
    DORADO (
        POD5.out.pod5,
        model
    )
    ch_versions = ch_versions.mix(DORADO.out.versions.first())

    emit:
    fastq    = DORADO.out.fastq      // channel: [ meta, fastq ]
    versions = ch_versions           // channel: [ path(versions.yml) ]
}
//
// Create depth and coverage statistics for each merged BAM file
//

include { SAMTOOLS_FLAGSTAT                    } from '../../../modules/nf-core/samtools/flagstat'
include { MOSDEPTH                             } from '../../../modules/nf-core/mosdepth'
include { SUMMARIZE_STATS                      } from '../../../modules/local/summarize_stats'

workflow BAM_STATS {
    take:
    ch_bam_bai     // channel (mandatory): [ val(meta), path(bam), path(bai) ]
    ch_fasta_fai   // channel (mandatory): [ val(meta), path(fasta), path(fai) ]
    ch_read_stats  // channel (mandatory): [ val(meta), path(json) ]
    bed_file       // file (ptional)

    main:
    ch_versions = Channel.empty()

    // get number of reads per sample before and after trimming:
    ch_read_stats.map { meta, json -> [ meta.sample, json ] }
        .groupTuple()
        .map { sample, jsons -> [ sample.toString(), WorkflowONTGeno.getNumberOfReads(jsons) ] }
        .set { ch_counts }

    // run SAMtools flagstat per merged bam file:
    SAMTOOLS_FLAGSTAT(ch_bam_bai)
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions.first())

    // run mosdepth per merged bam file:
    MOSDEPTH (
        ch_bam_bai.map { meta, bam, bai -> [ meta, bam, bai, bed_file ] },
        ch_fasta_fai.map{ meta, fasta, fai -> [ meta, fasta ] }
    )
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions.first())

    // Summarize all bam stats:
    SAMTOOLS_FLAGSTAT.out.flagstat
        .join(MOSDEPTH.out.summary_txt, failOnDuplicate:true, failOnMismatch:true)
        .join(MOSDEPTH.out.regions_txt, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, stat, depth, cov -> [ meta.sample, stat, depth, cov ] }
        .combine(ch_counts, by:0)
        .map { sample, stat, depth, cov, counts ->
            [ [id: 'summary'], stat, depth, cov, counts[0], counts[1] ]
          }
        .groupTuple()
        .set { ch_to_summary }

    SUMMARIZE_STATS(ch_to_summary)
    ch_versions = ch_versions.mix(SUMMARIZE_STATS.out.versions)

    emit:
    versions = ch_versions        // channel: [ versions.yml ]
}
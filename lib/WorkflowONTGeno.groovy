//
// This file holds several functions specific to the workflow/ontgeno.nf in the nf-core/ontgeno pipeline
//

import groovy.json.JsonSlurper

class WorkflowONTGeno {

    //
    // Function to parse FastP json files to get total number of reads before and after trimming per sample
    //
    public static Tuple getNumberOfReads(json_files) {
        def total_reads = 0
        def trimmed_reads = 0
        json_files.each { json_file ->
            def jsonSlurper = new JsonSlurper()
            def json = jsonSlurper.parse(json_file)
            total_reads += json.summary.before_filtering.total_reads
            trimmed_reads += json.summary.after_filtering.total_reads
        }
        return new Tuple(total_reads, trimmed_reads)
    }
}
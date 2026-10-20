nextflow.enable.dsl = 2
include { ANNOTATION_RELEASE } from './workflows/annotationRelease.nf'

workflow {
    def rows=new groovy.json.JsonSlurperClassic().parse(new File(params.release_test_requests))
    def references=Channel.fromList(rows).map { r ->
        tuple(r.aid,file(r.genome),file(r.masked),file(r.identity),file(r.repeats),file(r.repeat_status))
    }
    def identified=Channel.fromList(rows).map { r -> tuple(r.aid,file(r.identified)) }
    def functions=Channel.fromList(rows).filter { params.run_functional_annotation }.map { r ->
        tuple(r.aid,file(r.functional.gff),file(r.functional.table),file(r.functional.genes),file(r.functional.status))
    }
    def tracks=Channel.fromList(rows).flatMap { r ->
        r.tracks.collect { p -> tuple(r.aid,new File(p).name.split('\\.')[1],file(p)) }
    }
    def statuses=Channel.fromList(rows).map { r -> tuple(r.aid,file(r.track_status)) }
    ANNOTATION_RELEASE(references,identified,functions,tracks,statuses)
}

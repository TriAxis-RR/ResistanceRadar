nextflow.enable.dsl=2

include { FASTP } from './modules/fastp/main.nf'

params.samplesheet = "tests/samplesheets/erp015409_samples.csv"

workflow {
    reads_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            tuple(row.sample_id, file(row.read1), file(row.read2))
        }

    FASTP(reads_ch)
}

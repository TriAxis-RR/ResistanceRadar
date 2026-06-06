process FASTP {
    tag "$sample_id"

    cpus 4
    memory '8 GB'
    time '4h'
    publishDir "/rsstu/users/s/sleblan/MismatchRepair/Project/results/ERP015409/${sample_id}/nextflow_fastp", mode: 'copy'
    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id),
          path("${sample_id}_clean_R1.fastq.gz"),
          path("${sample_id}_clean_R2.fastq.gz"),
          path("${sample_id}_fastp.json"),
          path("${sample_id}_fastp.html")

    script:
    """
    fastp \
      -i ${read1} \
      -I ${read2} \
      -o ${sample_id}_clean_R1.fastq.gz \
      -O ${sample_id}_clean_R2.fastq.gz \
      --detect_adapter_for_pe \
      --qualified_quality_phred 20 \
      --unqualified_percent_limit 40 \
      --length_required 50 \
      --correction \
      --overlap_len_require 10 \
      --thread ${task.cpus} \
      --json ${sample_id}_fastp.json \
      --html ${sample_id}_fastp.html
    """
}

#! /usr/bin/env nextflow

ref = file(params.ref)
ref_fai = file(params.ref_fai)
ref_amb = file(params.ref_amb)
ref_ann = file(params.ref_ann)
ref_bwt = file(params.ref_bwt)
ref_pac = file(params.ref_pac)
ref_sa = file(params.ref_sa)
////ref_dict = file(params.ref_dict)
// fastq1Path = file(params.fastqr1)
// fastq2Path = file(params.fastqr2)
// sample_id = params.sample_id
threads = params.threads
outpath = params.outpath
code_ver = params.code_ver
iswes = params.wes
aws_region = params.aws_region
//jobdef = params.jobdef
reference_mills = file(params.ref_mills)
reference_mills_index = file(params.ref_mills_tbi)
reference_mills_idx = file(params.ref_mill_idx)
reference_dbsnp = file(params.ref_dbsnp)
reference_dbsnp_index = file(params.ref_dbsnp_idx)
//reference_bedfile = params.wes_bedfile


Channel
    .fromPath(params.input)
    .splitCsv(header:true)
    .map{ row-> tuple(row.sample_id, file(row.read1), file(row.read2)) }
    .into { samples1_ch; samples2_ch }

// process read_input_csv {
//input: 
    //  set sample_id, file(read1), file(read2) from samples1_ch 

// }

process alignment {
   //  tag "${sample_id}"
    stageInMode 'copy'
    cpus 128
    memory '245760 MB'	

    input:
   set sample_id, file(fastq1Path), file(fastq2Path) from samples1_ch
 file ref
//    file ref_dict
    file ref_fai
    file ref_amb
    file ref_ann
    file ref_bwt
    file ref_pac
    file ref_sa   
    // file fastq1Path
    // file fastq2Path    


    output:
    //  file  "bwa_${sample_id}_output.bam" into bwa_ch    
    //  file  "bwa_${sample_id}_output.bam.bai" into bwa_ch_bai    
    file "${sample_id}_sorted.bam" into (outputs_sorted_bam1, outputs_sorted_bam2,outputs_sorted_bam3)
    file "${sample_id}_sorted.bam.bai" into (outputs_indexed_bam1,outputs_indexed_bam2,outputs_indexed_bam3)
    val sample_id into (sample_id_ch_1, sample_id_ch_2,sample_id_ch_3,sample_id_ch_4,sample_id_ch_5,sample_id_ch_6,sample_id_ch_7,sample_id_ch_8,sample_id_ch_9,sample_id_ch_10,sample_id_ch_11,sample_id_ch_12,sample_id_ch_13)
    publishDir "${outpath}/fastq/${sample_id}/", mode: 'copy', overwrite: false
    script:
    """
    export SENTIEON_LICENSE=172.31.10.202:8990
    sentieon bwa mem -v 1 -K '${params.chunk_size}' -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:MGI" -M -t '${task.cpus}' '${ref}' '${fastq1Path}' '${fastq2Path}'| sentieon util sort -r '${ref}' -o '${sample_id}_sorted.bam' -t '${task.cpus}' --sam2bam -i -
    
    """
}
// process sentieon_util {
//    tag "${sample_id}"
//    cpus 10
//    input:
//        file local_alignment from bwa_ch 
//        file ref
 //       file local_alignment_bai from bwa_ch_bai 
//    output:
//        file "${sample_id}_sorted.bam" into (outputs_sorted_bam1, outputs_sorted_bam2,outputs_sorted_bam3)
//        file "${sample_id}_sorted.bam.bai" into (outputs_indexed_bam1,outputs_indexed_bam2,outputs_indexed_bam3)
//
//    script:
//    """
//    sentieon util sort -r '${ref}' -i local_alignment  -o '${sample_id}_sorted.bam' -t '${task.cpus}' 
//    """
// }

process alignment_metrics {
    tag "${sample_id}"
    stageInMode 'copy'
    cpus 128
    memory '245760 MB'	
    input:
    file aligned_bam from outputs_sorted_bam1
    file bam_index from outputs_indexed_bam1
    file ref
//    file ref_dict
    file ref_fai
    val sample_id from sample_id_ch_2
    output:
    file "${sample_id}.mq_metrics.txt" into output_meanqualitybycycle
    file "${sample_id}.qd_metrics.txt" into output_qualdistribution
    file "${sample_id}.gc_summary.txt" into output_gcsummary
    file "${sample_id}.gc_metrics.txt" into output_gcmetrics
    file "${sample_id}.aln_metrics.txt" into output_alignsummary
    file "${sample_id}.is_metrics.txt" into output_insertsize

    publishDir "${outpath}/metrics/${sample_id}/", mode: 'copy', overwrite: false

    script:
    """
    sentieon driver --input '${aligned_bam}' \
                    --reference '${ref}'  \
                    --thread_count '${task.cpus}' \
                    --algo MeanQualityByCycle '${sample_id}.mq_metrics.txt'  \
                    --algo QualDistribution '${sample_id}.qd_metrics.txt'  \
                    --algo GCBias --summary '${sample_id}.gc_summary.txt' '${sample_id}.gc_metrics.txt'  \
                    --algo AlignmentStat '${sample_id}.aln_metrics.txt'  \
                    --algo InsertSizeMetricAlgo '${sample_id}.is_metrics.txt'
    """
}
process plot_alignment_metrics_gc {
    tag "${sample_id}"
    errorStrategy 'finish'
    stageInMode 'copy'
    cpus 128
    memory '245760 MB'	
    input:
    file input_gcmetrics from output_gcmetrics
val sample_id from sample_id_ch_3

    output:
    file "${sample_id}.gc_metrics_plot.pdf"

    publishDir "${outpath}/metrics/${sample_id}/", mode: 'copy', overwrite: false

    script:
    """
    sentieon plot GCBias -o "${sample_id}.gc_metrics_plot.pdf" '${input_gcmetrics}'
    """
}

process plot_alignment_metrics_qd {
    tag "${sample_id}"
    stageInMode 'copy'
    cpus 128
    memory '245760 MB'	

    input:
    file input_qualdistribution from output_qualdistribution
val sample_id from sample_id_ch_4

    output:
    file "${sample_id}.qd_metrics_plot.pdf"

    publishDir "${outpath}/metrics/${sample_id}/", mode: 'copy', overwrite: false

    script:
    """
    sentieon plot QualDistribution -o "${sample_id}.qd_metrics_plot.pdf" '${input_qualdistribution}'
    """
}

process plot_alignment_metrics_mq {
    tag "${sample_id}"
    stageInMode 'copy'
    cpus 128
    memory '245760 MB'	

    input:
    file input_meanqualitybycycle from output_meanqualitybycycle
val sample_id from sample_id_ch_5

    output:
    file "${sample_id}.mq_metrics_plot.pdf"

    publishDir "${outpath}/metrics/${sample_id}/", mode: 'copy', overwrite: false

    script:
    """
    sentieon plot MeanQualityByCycle -o "${sample_id}.mq_metrics_plot.pdf" '${input_meanqualitybycycle}'
    """
}

process plot_alignment_metrics_isize {
    tag "${sample_id}"
    stageInMode 'copy'
    cpus 128
    memory '245760 MB'	

    input:
    file input_insertsize from output_insertsize
val sample_id from sample_id_ch_6

    output:
    file "${sample_id}.is_metrics_plot.pdf"

    publishDir "${outpath}/metrics/${sample_id}/", mode: 'copy', overwrite: false

    script:
    """
    sentieon plot InsertSizeMetricAlgo -o "${sample_id}.is_metrics_plot.pdf" '${input_insertsize}'
    """
}
process locus_collector {
    tag "${sample_id}"
    stageInMode 'copy'

    cpus 128
    memory '245760 MB'	

    input:
    file aligned_bam from outputs_sorted_bam2
    file bam_index from outputs_indexed_bam2
val sample_id from sample_id_ch_7

    output:
    file "${sample_id}.score.txt" into output_score_info
    file "${sample_id}.score.txt.idx" into output_score_info_index

    script:
    """
    sentieon driver --input '${aligned_bam}' \
                    --thread_count '${task.cpus}' \
                    --algo LocusCollector \
                    --fun score_info '${sample_id}.score.txt'
    """
}

process deduplication {
    tag "${sample_id}"
    stageInMode 'copy'
    cpus 128
    memory '245760 MB'	

    input:
    file aligned_bam from outputs_sorted_bam3
    file bam_index from outputs_indexed_bam3
    file score_info from output_score_info
    file score_info_index from output_score_info_index
val sample_id from sample_id_ch_8

    output:
    file "${sample_id}.deduped.bam" into (outputs_deduped_bam1,outputs_deduped_bam2,outputs_deduped_bam3)
    file "${sample_id}.deduped.bam.bai" into (outputs_deduped_indexed_bam1,outputs_deduped_indexed_bam2,outputs_deduped_indexed_bam3)
    file "${sample_id}.dedup_metrics.txt" into (outputs_deduped_metrics1,outputs_deduped_metrics2,outputs_deduped_metrics3)


    script:
    """
    sentieon driver --input '${aligned_bam}' \
                    --thread_count '${task.cpus}' \
                    --algo Dedup \
                    --metrics '${sample_id}.dedup_metrics.txt' \
                    --rmdup \
                    --score_info '${score_info}' \
                    '${sample_id}.deduped.bam'
    """
}
// process realignment {
//    tag "${sample_id}"

//  cpus 10

//    input:
//    file deduped_bam from outputs_deduped_bam
//    file bam_index from outputs_deduped_indexed_bam
//    file ref
////    file ref_dict
//    file ref_fai
//    file reference_mills
//	file reference_mills_index

//    output:
//    file "${sample_id}.bam" into outputs_bam_realignment
//    file "${sample_id}.bam.bai" into outputs_realigned_indexed_bam
//    file "code_ver"

 //   publishDir "${outpath}/bam/${sample_id}/", mode: 'copy', overwrite: false
//    script:
//    """
//    sentieon driver --input '${deduped_bam}' \
                        // --reference '${ref}' \
                        // --thread_count '${task.cpus}' \
                        // --algo Realigner \
                        // --known_sites  '${reference_mills}' \
                        // '${sample_id}.bam'

//    echo "${code_ver}" > code_ver
//    """
//}

process qualcal {
    tag "${sample_id}"
    stageInMode 'copy'
    cpus 128
    memory '245760 MB'	
    input:
    file deduped_bam from outputs_deduped_bam1
    file bam_index from outputs_deduped_indexed_bam1
    file ref
//    file ref_dict
    file ref_fai
    file reference_mills
    file reference_mills_index
    file reference_dbsnp
    file reference_dbsnp_index
    file reference_mills_idx
val sample_id from sample_id_ch_9

    output:
    file "${sample_id}.recal.table" into (outputs_recal_table1,outputs_recal_table2,outputs_recal_table3)

    script:
    """
    sentieon driver --input '${deduped_bam}' \
                        --reference '${ref}' \
                        --thread_count '${task.cpus}' \
                        --algo QualCal \
                        --known_sites '${reference_dbsnp}' \
                        --known_sites '${reference_mills}' \
                        '${sample_id}.recal.table'
    """
}

process qualcalpost {
    tag "${sample_id}"
    stageInMode 'copy'
    cpus 128
    memory '245760 MB'	

    input:
    file deduped_bam from outputs_deduped_bam2
    file bam_index from outputs_deduped_indexed_bam2
    file qualcal from outputs_recal_table1
    file ref
//    file ref_dict
    file ref_fai
    file reference_mills
    file reference_mills_index
    file reference_dbsnp
    file reference_dbsnp_index
val sample_id from sample_id_ch_10

    file reference_mills_idx
    output:
    file "${sample_id}.recal.table.post" into outputs_recal_table_post

    script:
    """
    sentieon driver --input '${deduped_bam}' \
                        --qual_cal '${qualcal}' \
                        --reference '${ref}' \
                        --thread_count '${task.cpus}' \
                        --algo QualCal \
                        --known_sites '${reference_dbsnp}' \
                        --known_sites '${reference_mills}'\
                        '${sample_id}.recal.table.post'
    """
}

process applyrecal {
    tag "${sample_id}"
    stageInMode 'copy'
    cpus 128
    memory '245760 MB'	

    input:
    file qualcal from outputs_recal_table2
    file recal_table from outputs_recal_table_post
val sample_id from sample_id_ch_11

    output:
    file "${sample_id}.recal.csv" into outputs_bqsr

    script:
    """
    sentieon driver --thread_count '${task.cpus}' \
                        --algo QualCal \
                        --after '${recal_table}' \
                        --before '${qualcal}' \
                        --plot '${sample_id}.recal.csv'
    """
}

process plotbqsr {
    tag "${sample_id}"
    stageInMode 'copy'
    cpus 128
    memory '245760 MB'	

    input:
    file recal_table from outputs_bqsr
val sample_id from sample_id_ch_12

    output:
    file "${sample_id}.recal.pdf"

    publishDir "${outpath}/metrics/${sample_id}/", mode: 'copy', overwrite: false

    script:
    """
    sentieon plot QualCal -o '${sample_id}.recal.pdf' '${recal_table}'
    """
}

process haplotypecaller {
    tag "${sample_id}"
    stageInMode 'copy'
    cpus 128
    memory '245760 MB'	

    input:
    file deduped_bam from outputs_deduped_bam3
    file bam_index from outputs_deduped_indexed_bam3
    file recal_table from outputs_recal_table3
    file ref
//    file ref_dict
    file ref_fai
    file reference_dbsnp
    file reference_dbsnp_index
val sample_id from sample_id_ch_13

    output:
    file "${sample_id}.g.vcf.gz" into outputs_gvcf
    file "${sample_id}.g.vcf.gz.tbi" into outputs_gvcf_index
    

    publishDir "${outpath}/gvcf/${sample_id}/", mode: 'copy', overwrite: false

   
  //  ARG_INTERVAL=""
//    if [[ ${iswes} == "true" ]]
//    then
//        aws s3 cp ${reference_bedfile} wes.bed --region ${aws_region}
//        ARG_INTERVAL="--interval wes.bed"
//    fi
    script:
    """   
    sentieon driver  --input '${deduped_bam}' \
                    --qual_cal '${recal_table}' \
                    --reference '${ref}' \
                    --thread_count '${task.cpus}' \
                    --algo Haplotyper \
                    --call_conf 30 \
                    --dbsnp '${reference_dbsnp}' \
                    --emit_conf 30 \
                    --emit_mode Gvcf \
                    '${sample_id}.g.vcf.gz'

    """

}


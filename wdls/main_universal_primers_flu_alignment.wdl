version 1.0

import "./tasks_alignment.wdl" as tasks

workflow universal_primers_flu_alignment {

    input{
        String sample_name
        String project_name
        String out_dir # base project path (don't include terra outputs)
        # e.g. gs://..../flu/flu_0000_miseq/realingment/

        String subtype # this will be the subdirectory name

        File fastq1
        File fastq2

        # reference sequences
        File reference 

        # python script
        File summary_python_script

    }

    # private declaration
    # docker containers
    String fastqc_docker = "staphb/fastqc:0.12.1"
    String fastp_docker = "staphb/fastp:0.23.2"
    String ivar_docker = "staphb/ivar:1.4.4-aligners" # has bwa and samtools
    String ubuntu_docker = 'ubuntu:jammy-20240627.1'
    String utility_docker = 'theiagen/utility:1.0'
    String python_docker = "mchether/py3-bio:v1"

    # workflow version
    String workflow_version = "v0_1_0"

    # transfer_path
    String transfer_path = "~{out_dir}/~{workflow_version}/~{subtype}/"
    # e.g. gs://...../flu/flu_0000_miseq/realingment/v0.0.0/{subtype}/

    call tasks.fastqc as fastqc_raw{
        input:
            sample_name = sample_name,
            fastq1 = fastq1,
            fastq2 = fastq2,
            fastq_type = "raw",
            docker = fastqc_docker
    }
    
    call tasks.fastp {
        input:
            sample_name = sample_name,
            fastq1 = fastq1,
            fastq2 = fastq2,
            docker = fastp_docker
    }
    call tasks.fastqc as fastqc_clean{
        input:
            sample_name = sample_name,
            fastq1 = fastp.fastq_1_cleaned,
            fastq2 = fastp.fastq_2_cleaned,
            fastq_type = "cleaned",
            docker = fastqc_docker
    }

    call tasks.align_bwa {
        input:
            sample_name = sample_name,
            fastq_1_cleaned = fastp.fastq_1_cleaned,
            fastq_2_cleaned = fastp.fastq_2_cleaned,
            reference = reference,
            docker = ivar_docker
    }

    call tasks.call_consensus_ivar {
        input:
            sample_name = sample_name,
            bam = align_bwa.bam,
            reference = reference,
            docker = ivar_docker
    }

    call tasks.calculate_alignment_metrics {
        input:
            sample_name = sample_name,
            bam = align_bwa.bam,
            reference = reference,
            docker = ivar_docker

    }

    call tasks.summarize_results {
        input:
            sample_name = sample_name,
            project_name = project_name,
            consensus_fasta = call_consensus_ivar.consensus_fasta,
            reference = reference,
            fastqc_clean_summary = fastqc_clean.summary_metrics,
            fastqc_raw_summary = fastqc_raw.summary_metrics,
            samtools_coverage = calculate_alignment_metrics.coverage,
            summary_python_script = summary_python_script,
            docker = python_docker
    }

    call tasks.version_capture {
        input:
            sample_name = sample_name,
            project_name = project_name,
            workflow_version = workflow_version,
            fastqc_version = fastqc_raw.fastqc_version,
            fastp_version = fastp.fastp_version,
            bwa_version = align_bwa.bwa_version,
            ivar_version = call_consensus_ivar.ivar_version,
            samtools_version = align_bwa.samtools_version,
            docker = ubuntu_docker
    }

    call tasks.transfer {
        input:
            transfer_path = transfer_path,
            fastqc1_data_raw = fastqc_raw.fastqc1_data,
            fastqc2_data_raw = fastqc_raw.fastqc2_data,
            fastqc1_html_raw = fastqc_raw.fastqc1_html,
            fastqc2_html_raw = fastqc_raw.fastqc2_html,
            fastqc1_data_clean = fastqc_clean.fastqc1_data,
            fastqc2_data_clean = fastqc_clean.fastqc2_data,
            fastqc1_html_clean = fastqc_clean.fastqc1_html,
            fastqc2_html_clean = fastqc_clean.fastqc2_html,
            fastp_fastq_1_cleaned = fastp.fastq_1_cleaned,
            fastp_fastq_2_cleaned = fastp.fastq_2_cleaned,
            fastp_html = fastp.fastp_html,
            fastp_json = fastp.fastp_json,
            bam = align_bwa.bam,
            bai = align_bwa.bai,
            consensus_fasta = call_consensus_ivar.consensus_fasta,
            coverage = calculate_alignment_metrics.coverage,
            stats = calculate_alignment_metrics.stats,
            depth = calculate_alignment_metrics.depth,
            summary_results = summarize_results.summary_results,
            version_file = version_capture.version_file,
            docker = utility_docker



    }

    output {

        File fastqc1_data_raw = fastqc_raw.fastqc1_data
        File fastqc2_data_raw = fastqc_raw.fastqc2_data
        File fastqc1_html_raw = fastqc_raw.fastqc1_html
        File fastqc2_html_raw = fastqc_raw.fastqc2_html
        String fastqc_summary_metrics_raw = fastqc_raw.summary_metrics

        File fastqc1_data_clean = fastqc_clean.fastqc1_data
        File fastqc2_data_clean = fastqc_clean.fastqc2_data
        File fastqc1_html_clean = fastqc_clean.fastqc1_html
        File fastqc2_html_clean = fastqc_clean.fastqc2_html
        String fastqc_summary_metrics_clean = fastqc_clean.summary_metrics

        File fastp_fastq_1_cleaned = fastp.fastq_1_cleaned
        File fastp_fastq_2_cleaned = fastp.fastq_2_cleaned
        File fastp_html = fastp.fastp_html
        File fastp_json = fastp.fastp_json

        File bam = align_bwa.bam
        File bai = align_bwa.bai

        File consensus_fasta = call_consensus_ivar.consensus_fasta

        File coverage = calculate_alignment_metrics.coverage
        File stats = calculate_alignment_metrics.stats
        File depth = calculate_alignment_metrics.depth

        File summary_results = summarize_results.summary_results
        File version_file = version_capture.version_file
    }
}



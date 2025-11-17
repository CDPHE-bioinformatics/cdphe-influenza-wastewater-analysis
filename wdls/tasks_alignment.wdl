version 1.0

task fastqc {
    input {
        String sample_name
        File fastq1
        File fastq2
        String fastq_type
        String docker
    }

    String fastq1_name = basename(fastq1, ".fastq.gz")
    String fastq2_name = basename(fastq2, ".fastq.gz")
    String fastq1_string = "~{fastq1_name}_fastqc"
    String fastq2_string = "~{fastq2_name}_fastqc"
    String fastq1_data = "~{fastq1_string}_data.txt"
    String fastq2_data = "~{fastq2_string}_data.txt"
    String fastq1_html = "~{fastq1_string}_report.html"
    String fastq2_html = "~{fastq2_string}_report.html"
    String summary_metrics_fn = "~{sample_name}_~{fastq_type}_summary_metrics.csv"

    command <<<
        echo "Running FastQC on ~{fastq_type} reads for sample ~{sample_name}"

        fastqc --outdir $PWD --extract --delete ~{fastq1} ~{fastq2}
        fastqc --version | awk '{print $2}' | tee VERSION  
        mv "~{fastq1_string}/fastqc_data.txt" ~{fastq1_data}
        mv "~{fastq1_string}/fastqc_report.html" ~{fastq1_html}
        mv "~{fastq2_string}/fastqc_data.txt" ~{fastq2_data}
        mv "~{fastq2_string}/fastqc_report.html" ~{fastq2_html}

        # Summarize output to simpler csv file
        summarize_fastqc () {
            total_seqs=$(grep "Total Sequences" $1 | cut -f 2)
            flagged_reads=$(grep "Sequences flagged as poor quality" $1 | cut -f 2)
            sequence_length=$(grep "Sequence length" $1 | cut -f 2)
            echo $total_seqs,$flagged_reads,$sequence_length
        }

        fastq1_data="~{fastq1_data}"
        fastq2_data="~{fastq2_data}"
        r1_info=$(summarize_fastqc ${fastq1_data})
        r2_info=$(summarize_fastqc ${fastq2_data})
        echo "sample_name,r1_total_reads,r1_flagged_reads_as_poor_quality,r1_read_len,r2_total_reads,r2_flagged_reads_as_poor_quality,r2_read_len" >> ~{summary_metrics_fn} 
        echo "~{sample_name},${r1_info},${r2_info}" >> ~{summary_metrics_fn}
    >>>

    output{
        File fastqc1_data = "~{fastq1_data}"
        File fastqc1_html = "~{fastq1_html}"
        File fastqc2_data = "~{fastq2_data}"
        File fastqc2_html = "~{fastq2_html}"
        File summary_metrics = "~{summary_metrics_fn}"
        String fastqc_version = read_string('VERSION')
    }

    runtime {
        cpu: 2
        memory: "2G"
        disks: "local-disk 2 HDD"
        docker: docker
    }

}

task fastp {
    input{
        String sample_name
        File fastq1
        File fastq2
        String docker
    }
    
    String fastq1_name = basename(fastq1, ".fastq.gz")
    String fastq2_name = basename(fastq2, ".fastq.gz")
    String cleaned_1_name = "~{fastq1_name}_clean.fastq.gz"
    String cleaned_2_name = "~{fastq2_name}_clean.fastq.gz"
    String unpaired_1_name = "~{fastq1_name}_unpaired.fastq.gz"
    String unpaired_2_name = "~{fastq2_name}_unpaired.fastq.gz"
    String fastp_html_name = "~{sample_name}_fastp.html"
    String fastp_json_name = "~{sample_name}_fastp.json"

    command <<<
        fastp --version 2>&1 | awk '{print $2}' | tee VERSION
        fastp \
            --in1 ~{fastq1} --in2 ~{fastq2} \
            --out1 ~{cleaned_1_name} --out2 ~{cleaned_2_name}  \
            --unpaired1 ~{unpaired_1_name} --unpaired2 ~{unpaired_2_name} \
            --cut_tail \
            --cut_tail_window_size 4 \
            --cut_tail_mean_quality 30 \
            --length_required 70 \
            --detect_adapter_for_pe \
            --trim_poly_g \
            --html ~{fastp_html_name} \
            --json ~{fastp_json_name}
    >>>

    output {
        File fastq_1_cleaned = cleaned_1_name
        File fastq_2_cleaned = cleaned_2_name
        File fastq_1_unpaired = unpaired_1_name
        File fastq_2_unpaired = unpaired_2_name
        File fastp_html = fastp_html_name
        File fastp_json = fastp_json_name
        String fastp_version = read_string('VERSION')

    }

    runtime {
        docker: docker
        memory: "8 GB"
        cpu:    4
        disks: "local-disk 50 SSD"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
}

task align_bwa {
    input {
        String sample_name
        File fastq_1_cleaned
        File fastq_2_cleaned
        File reference
        String docker

    }

    String reference_name = basename(reference, ".fasta")
    String sam_fn = "~{sample_name}.sam"
    String bam_fn = "~{sample_name}.bam"
    String bai_fn = "~{sample_name}.bam.bai"

    command <<<

        samtools --version-only | tee SAMTOOLS_VERSION
        bwa 2>&1 | awk '/Version/{print $2}' | tee BWA_VERSION

        bwa index -p ~{reference_name} -a is ~{reference}
        bwa mem -t 4 ~{reference_name} ~{fastq_1_cleaned} ~{fastq_2_cleaned} > ~{sam_fn}
        samtools view -b -@ 4 ~{sam_fn} | samtools sort -m 2G -@ 4 -o ~{bam_fn}
        samtools index ~{bam_fn} 
    >>>

    output {
        File bam = bam_fn
        File bai = bai_fn
        String bwa_version = read_string('BWA_VERSION')
        String samtools_version = read_string('SAMTOOLS_VERSION')

    }

    runtime {
        cpu: 2
        memory: "2G"
        disks: "local-disk 2 HDD"
        docker: docker
    }
}

task call_consensus_ivar {

    input {
        String sample_name
        File bam
        File reference
        String docker
    }

    String pileup_fn = "~{sample_name}_pileup.txt"
    String consensus_fn_prefix = "~{sample_name}_consensus"
    String consensus_fn = consensus_fn_prefix + ".fa"

    command <<<
        ivar version | awk '/version/ {print $3}' | tee VERSION
        samtools mpileup -A -aa -d 600000 -B -Q 20 -q 20 -f ~{reference} ~{bam} | \
        ivar consensus -p ~{consensus_fn_prefix} -q 20 -t 0.6 -m 10

        # not renaming fasta header becuase I want it to be compatible with the H5 workflow
    >>>

    output {
        File consensus_fasta = "~{consensus_fn}"
        String ivar_version = read_string('VERSION')
    }

    runtime {
        cpu: 2
        memory: "2G"
        disks: "local-disk 2 HDD"
        docker: docker
    }

}

task calculate_alignment_metrics {
    input {
        String sample_name
        File bam
        File reference
        String docker
    }

    String coverage_fn = "~{sample_name}_coverage.txt"
    String stats_fn = "~{sample_name}_stats.txt"
    String depth_fn = "~{sample_name}_depth.txt"

    command <<<
        samtools depth -a -o ~{depth_fn} ~{bam}
        samtools coverage -o ~{coverage_fn} ~{bam}
        samtools stats ~{bam} > ~{stats_fn}
    >>>

    output {
        File coverage = "~{coverage_fn}"
        File stats = "~{stats_fn}"
        File depth = "~{depth_fn}"
    }

    runtime {
        cpu: 2
        memory: "4 GB"
        docker: docker
    }
}

task summarize_results {
    input {
        String sample_name
        String project_name
        File consensus_fasta
        File reference
        File fastqc_clean_summary
        File fastqc_raw_summary
        File samtools_coverage
        File summary_python_script
        String docker
    }

    String summary_results_fn = "~{sample_name}_results.csv"

    command <<<
        python ~{summary_python_script} \
            --sample_name ~{sample_name} \
            --project_name ~{project_name} \
            --consensus_fasta ~{consensus_fasta} \
            --reference ~{reference} \
            --fastqc_clean_summary ~{fastqc_clean_summary} \
            --fastqc_raw_summary ~{fastqc_raw_summary} \
            --samtools_coverage ~{samtools_coverage} \
            --out_fn ~{summary_results_fn}
    >>>

    output {
        File summary_results = "~{summary_results_fn}"
    }

    runtime {
        cpu: 1
        memory: "2 GB"
        docker: docker
    }
}

task version_capture{

    input {
        String sample_name
        String project_name
        String workflow_version
        String fastqc_version
        String fastp_version
        String bwa_version
        String ivar_version
        String samtools_version
        String docker
    }

    String version_fn = "~{sample_name}_versions.txt"

    command <<<
        echo "Sample Name: ~{sample_name}" > ~{version_fn}
        echo "Project Name: ~{project_name}" >> ~{version_fn}
        echo "Workflow Version: ~{workflow_version}" >> ~{version_fn}
        echo "FastQC Version: ~{fastqc_version}" >> ~{version_fn}
        echo "Fastp Version: ~{fastp_version}" >> ~{version_fn}
        echo "BWA Version: ~{bwa_version}" >> ~{version_fn}
        echo "iVar Version: ~{ivar_version}" >> ~{version_fn}
        echo "Samtools Version: ~{samtools_version}" >> ~{version_fn}

    >>>

    output {  
        File version_file = "~{version_fn}"
    }

    runtime {
        docker: docker
        memory: "2 GB"
        cpu:    1
    }
}

task transfer {

    input {
        String transfer_path

        File fastqc1_data_raw
        File fastqc2_data_raw
        File fastqc1_html_raw
        File fastqc2_html_raw

        File fastqc1_data_clean
        File fastqc2_data_clean
        File fastqc1_html_clean
        File fastqc2_html_clean

        File fastp_fastq_1_cleaned
        File fastp_fastq_2_cleaned
        File fastp_html
        File fastp_json

        File bam
        File bai

        File consensus_fasta

        File coverage
        File stats
        File depth

        File summary_results
        File version_file

        String docker

    }

    command <<<
        gsutil -m cp -r ~{fastqc1_data_raw} ~{transfer_path}/fastqc_raw/
        gsutil -m cp -r ~{fastqc2_data_raw} ~{transfer_path}/fastqc_raw/
        gsutil -m cp -r ~{fastqc1_html_raw} ~{transfer_path}/fastqc_raw/
        gsutil -m cp -r ~{fastqc2_html_raw} ~{transfer_path}/fastqc_raw/
    
        gsutil -m cp -r ~{fastqc1_data_clean} ~{transfer_path}/fastqc_clean/
        gsutil -m cp -r ~{fastqc2_data_clean} ~{transfer_path}/fastqc_clean/
        gsutil -m cp -r ~{fastqc1_html_clean} ~{transfer_path}/fastqc_clean/
        gsutil -m cp -r ~{fastqc2_html_clean} ~{transfer_path}/fastqc_clean/

        gsutil -m cp -r ~{fastp_fastq_1_cleaned} ~{transfer_path}/fastp/
        gsutil -m cp -r ~{fastp_fastq_2_cleaned} ~{transfer_path}/fastp/
        gsutil -m cp -r ~{fastp_html} ~{transfer_path}/fastp/
        gsutil -m cp -r ~{fastp_json} ~{transfer_path}/fastp/

        gsutil -m cp -r ~{bam} ~{transfer_path}/alignment/
        gsutil -m cp -r ~{bai} ~{transfer_path}/alignment/

        gsutil -m cp -r ~{consensus_fasta} ~{transfer_path}/consensus/

        gsutil -m cp -r ~{coverage} ~{transfer_path}/samtools/
        gsutil -m cp -r ~{stats} ~{transfer_path}/samtools/
        gsutil -m cp -r ~{depth} ~{transfer_path}/samtools/

        gsutil -m cp -r ~{summary_results} ~{transfer_path}/summary_results/
        gsutil -m cp -r ~{version_file} ~{transfer_path}/version_capture/

    >>>

    runtime {
        docker: docker
        memory: "2 GB"
        cpu:    1
    }
}


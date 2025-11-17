version 1.0

task freyja{
    input {
        String sample_name
        File bam

        # reference sequences
        File reference
        File gff
        File barcodes
        String docker
    }
        String variants_fn = "~{sample_name}_variants.tsv"
        String depths_fn = "~{sample_name}_depths.depth"
        String demix_fn = "demix/~{sample_name}_demix.tsv"
        String aggregate_fn = "~{sample_name}_aggregate.tsv"

    command {

        # grab freyja version
        freyja --version | cut -d',' -f2 | tee VERSION

        # run variants to get variants file and depth file
        mkdir demix # have to make a demix directory to call aggregrate since input is directory

        freyja variants ~{bam} --ref ~{reference} --annot ~{gff} --variants ~{variants_fn} --depths ~{depths_fn}
        freyja demix ~{variants_fn} ~{depths_fn} --barcodes ~{barcodes} --output ~{demix_fn}
        freyja aggregate demix/ --output ~{aggregate_fn}
        
    }

    output{
        File variants_file = variants_fn
        File depths_file = depths_fn
        File demix_file = demix_fn
        File aggregate_file = aggregate_fn
        String freyja_version = read_string("VERSION")
    }

    runtime {
        docker: docker
        memory: "16 GB"
        cpu: 4
    }
}


task format_results {
    input {
        File format_aggregate_py_script
        File aggreagate_file
        String sample_name
        String docker
    }

    String out_fn = "~{sample_name}_formatted_aggregate.csv"

    command {
        python ~{format_aggregate_py_script} \
        --aggregate ~{aggreagate_file} \
        --sample_name ~{sample_name} \
        --out_fn ~{out_fn}   
    }

    output {
        File formatted_aggregate = "~{out_fn}"
    }

    runtime {
        docker: docker
        memory: "16 GB"
        cpu:    4
        disks: "local-disk 100 SSD"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
}

task version_capture {
    input {
        String sample_name
        String project_name
        String freyja_version
        String workflow_version
        String barcodes_version
        String docker
    }

    command {
        echo "Sample Name: ~{sample_name}" > ~{sample_name}_version.txt
        echo "Project Name: ~{project_name}" >> ~{sample_name}_version.txt
        echo "Workflow Version: ~{workflow_version}" >> ~{sample_name}_version.txt
        echo "Freyja Version: ~{freyja_version}" >> ~{sample_name}_version.txt
        echo "Freyja Barcodes Version: ~{barcodes_version}" >> ~{sample_name}_version.txt
        
    }

    output {
        File version_file = "~{sample_name}_version.txt"
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

        File freyja_variants
        File freyja_depths
        File freyja_demix
        File freyja_aggregate
        File formatted_aggregate
        File version_file

        String docker

    }

    command {
        gsutil -m cp -r ~{freyja_variants} ~{transfer_path}/variants/
        gsutil -m cp -r ~{freyja_depths} ~{transfer_path}/depths/
        gsutil -m cp -r ~{freyja_demix} ~{transfer_path}/demix/
        gsutil -m cp -r ~{freyja_aggregate} ~{transfer_path}/aggregate/
        gsutil -m cp -r ~{formatted_aggregate} ~{transfer_path}/formatted_aggregate/
        gsutil -m cp -r ~{version_file} ~{transfer_path}/version_capture/

    }

    output {
        String transfer_complete = "Transfer to ~{transfer_path} completed."
    }

    runtime {
        docker: docker
        memory: "2 GB"
        cpu:    1
    }
}
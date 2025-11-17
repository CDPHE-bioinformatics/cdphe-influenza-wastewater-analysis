version 1.0
import "./tasks_variant_calling.wdl" as tasks

workflow freyja_other_pathogens {

    input {
        String sample_name
        String project_name
        String out_dir # base project path (don't include terra outputs)
        # e.g. gs://..../flu/flu_0000_miseq/wwt_variant_calling/

        String subtype # this will be the subdirectory 
        File bam
        File reference
        File gff
        File barcodes
        String barcodes_version

        File format_aggregate_py_script

        
    }

    # private declarations
    # docker container
    String freyja_docker = "staphb/freyja:2.0.1-11_02_2025-00-37-2025-11-03"
    String python_docker = "mchether/py3-bio:v2"
    String utility_docker = 'theiagen/utility:1.0'
    String ubuntu_docker = 'ubuntu:jammy-20240627.1'

    # workflow version
    String workflow_version = "v0_1_0"

    # transfer path
    String transfer_path = "~{out_dir}/~{workflow_version}/~{subtype}"
    # e.g. gs://...../flu/flu_0000_miseq/wwt_variant_calling/v0.0.0/{subtype}/

    
    call tasks.freyja {
        input:
            sample_name = sample_name,
            bam = bam,
            reference = reference,
            gff = gff,
            barcodes = barcodes,
            docker = freyja_docker


    }

    call tasks.format_results {
        input:
            format_aggregate_py_script = format_aggregate_py_script,
            aggreagate_file = freyja.aggregate_file,
            sample_name = sample_name,
            docker = python_docker

    }

    call tasks.version_capture {
        input:
            sample_name = sample_name,
            project_name = project_name,
            freyja_version = freyja.freyja_version,
            workflow_version = workflow_version,
            barcodes_version = barcodes_version,
            docker = ubuntu_docker


    }

    call tasks.transfer {
        input:
            transfer_path = transfer_path,
            freyja_variants = freyja.variants_file,
            freyja_depths = freyja.depths_file,
            freyja_demix = freyja.demix_file,
            freyja_aggregate = freyja.aggregate_file,
            formatted_aggregate = format_results.formatted_aggregate,
            version_file = version_capture.version_file,
            docker = utility_docker

    }

    output {
        File freyja_variants = freyja.variants_file
        File freyja_depths = freyja.depths_file
        File freyja_demix = freyja.demix_file
        File freyja_aggregate = freyja.aggregate_file
        File formatted_aggregate = format_results.formatted_aggregate
        File version_file = version_capture.version_file

    }

}


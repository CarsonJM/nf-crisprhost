#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CarsonJM/nf-phist
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/CarsonJM/nf-phist
----------------------------------------------------------------------------------------
    Overview:
        1. Split input files into chunks of X genomes (Nextflow)
        2. Download genome chunks (process - aria2c)
        3. Run downloaded chunks through PHIST (process - phist)
        4. Delete downloaded chunks that have been run through phist (process - rm)
        5. Combine PHIST outputs from each chunk into one file (process)
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DEFINE FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
process ARIA2C {
    tag "${meta.id}"
    label "process_medium"
    storeDir "tmp/aria2c/${meta.id}"

    conda "envs/aria2c.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/aria2:1.36.0' :
        'biocontainers/aria2:1.36.0' }"

    input:
    tuple val(meta), val(urls)

    output:
    tuple val(meta), path("host_fastas/")       , emit: host_fastas
    tuple val(meta), path("download_complete")  , emit: download_complete

    script:
    def download_list   = urls.collect { url -> url.toString() }.join(',')
    """
    # create an input file for aria2c
    IFS=',' read -r -a download_array <<< "${download_list}"
    printf '%s\\n' "\${download_array[@]}" > aria2_file.tsv

    # download fasta files with aria2c
    aria2c \\
        --input=aria2_file.tsv \\
        --dir=host_fastas \\
        --max-concurrent-downloads=${task.cpus}

    touch download_complete
    """
}

process PHIST {
    tag "${meta.id}"
    label "process_high"
    storeDir "tmp/phist/${meta.id}"
    containerOptions "${ workflow.containerEngine == 'singularity' ?
        '-B ' + workflow.launchDir :
        '' }"

    conda "envs/phist.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b7/b72245719494da16ebe64f0f73e8f29f9880bb89dce1065a5d30c72b82f7cf68/data' :
        'community.wave.seqera.io/library/kmer-db_python:2fcd54c55e4e0870' }"

    input:
    tuple val(meta), val(host_fastas)
    path(virus_fasta)

    output:
    tuple val(meta), path("${meta.id}.phist_table.csv") , emit: phist_tables

    script:
    """
    # run phist on virus fasta and host fasta chunk
    phist.py \\
        ${virus_fasta} \\
        ${host_fastas} \\
        ${meta.id}.phist_table.csv \\
        ${meta.id}.phist_preds.csv \\
        -t ${task.cpus}
    """
}

process RM_GENOMES {
    tag "${meta.id}"
    label "process_single"
    storeDir "tmp/rm_genomes/${meta.id}"

    input:
    tuple val(meta) , val(host_fastas)
    tuple val(meta2), val(phist_tables)

    output:
    path("rm_complete")

    script:
    """
    # delete downloaded fastas that have been run through phist
    rm -rf ${host_fastas}/*

    touch rm_complete
    """
}

process COMBINE_PHIST {
    label "process_single"
    storeDir "."

    input:
    path(phist_tables)

    output:
    path("${params.output}")    , emit: final_output

    script:
    """
    # iterate over phist tables
    for table in ${phist_tables[0]}; do
       head -n 2 \${table} > ${params.output}
    done

    for table in ${phist_tables}; do
        tail -n +3 \${table} >> ${params.output}
    done
    """
}

// Run entry workflow
workflow {
    main:
    // Check if output file already exists
    def output_file = file("${params.output}")
    if (!output_file.exists()) {

        // 1. Split input files into chunks of X genomes (Nextflow)
        ch_virus_fasta = channel.fromPath(params.virus_fasta).collect()
        ch_host_fastas = channel.fromPath(params.host_file)
            .splitCsv(header: false, strip: true)
            .flatten()
            .collate(params.chunk_size)
            .toList()
            .flatMap{ file -> file.withIndex() }
            .map { file, index ->
                [ [ id: 'chunk_' + index ], file ]
            }

        // 2. Download genome chunks (process - aria2c)
        ARIA2C(
            ch_host_fastas
        )

        // 3. Run downloaded chunks through PHIST (process - phist)
        PHIST(
            ARIA2C.out.host_fastas,
            ch_virus_fasta
        )

        // 4. Delete downloaded chunks that have been run through phist (process - rm)
        RM_GENOMES(
            ARIA2C.out.host_fastas,
            PHIST.out.phist_tables
        )

        // 5. Combine PHIST outputs from each chunk into one file (process)
        COMBINE_PHIST(
            PHIST.out.phist_tables.map { _meta, tables -> [ tables ] }.collect()
        )
    } else {
        println "Output file [${params.output}] already exists! Skipping nf-phist."
    }

    // Delete intermediate and Nextflow-specific files
    workflow.onComplete {
        if (output_file.exists()) {
            def work_dir = new File("./work/")
            def tmp_dir = new File("./tmp/")
            def nextflow_dir = new File("./.nextflow/")
            def launch_dir = new File(".")

            work_dir.deleteDir()
            tmp_dir.deleteDir()
            nextflow_dir.deleteDir()
            launch_dir.eachFileRecurse { file ->
                if (file.name ==~ /\.nextflow\.log.*/) {
                    file.delete()
                }
            }
        }
    }
}


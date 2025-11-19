#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CarsonJM/nf-crisprhost
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/CarsonJM/nf-crisprhost
----------------------------------------------------------------------------------------
    Overview:
        1. Download CRISPR spacer file (process - Nextflow)
        2. Split CRISPR spacer file into chunks (process - SEQKIT_SPLIT2)
        3. Create spacerextractor DB from virus fasta (process - SPACEREXTRACTOR_CREATETARGETDB)
        4. Align CRISPR chunks to spacerextractor DB (process - SPACEREXTRACTOR_MAPTOTARGET)
        5. Combine all mapping results into single output file (process - SPACEREXTRACTOR_COMBINEMAPS)
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DEFINE FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
process SEQKIT_SPLIT2 {
    label 'process_super_high'
    storeDir "tmp/seqkit_split2"

    conda "envs/seqkit.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0' :
        'biocontainers/seqkit:2.9.0--h9ee0642_0' }"

    input:
    path(fasta)

    output:
    path("spacer_split/*")  , emit: split_fastas

    script:
    """
    seqkit \\
        split2 \\
            ${fasta} \\
            --threads ${task.cpus} \\
            --by-size ${params.chunk_size} \\
            --out-dir spacer_split
    """
}

process SPACEREXTRACTOR_CREATETARGETDB {
    label 'process_high'
    storeDir "tmp/spacerextractor/createtargetdb"

    conda "envs/spacerextractor.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/29/2900c330cd8dd25094b4cd86e4a32a576ddb340f412ac17f6715ac4136cf495c/data' :
        'biocontainers/spacerextractor:0.9.7--pyhdfd78af_0' }"

    input:
    path(virus_fasta)

    output:
    path("virus_targets_db/")   , emit: db

    script:
    def is_compressed = virus_fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = virus_fasta.getName().replace(".gz", "")
    """
    # if gzipped, decompress virus fasta
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${virus_fasta} > ${fasta_name}
    fi
    
    # create spacerextractor target db
    spacerextractor \\
        create_target_db \\
            -i ${fasta_name} \\
            -d virus_targets_db \\
            -t ${task.cpus} \\
            --replace_spaces
    """
}

process SPACEREXTRACTOR_MAPTOTARGET {
    tag "${meta.id}"
    label 'process_medium'
    storeDir "tmp/spacerextractor/maptotarget/${meta.id}"

    conda "envs/spacerextractor.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/29/2900c330cd8dd25094b4cd86e4a32a576ddb340f412ac17f6715ac4136cf495c/data' :
        'biocontainers/spacerextractor:0.9.7--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(spacer_fasta)
    path(target_db)

    output:
    tuple val(meta), path("${meta.id}_all_hits.tsv")    , emit: mapped_results

    script:
    """
    # create new workdir for pre-emption
    nxf_workdir=\$(pwd)
    mkdir -p ${workflow.workDir.toAbsolutePath()}/spacerextractor_maptotarget/${meta.id}/
    cd ${workflow.workDir.toAbsolutePath()}/spacerextractor_maptotarget/${meta.id}/

    # run spacerextractor mapping
    SE_map_get_hits.py \\
        map_to_target \\
            -i \${nxf_workdir}/${spacer_fasta} \\
            -d \${nxf_workdir}/${target_db} \\
            -o ${meta.id}_map_results \\
            -t ${task.cpus}

    # move tsv file to cwd
    mv ${meta.id}_map_results/${meta.id}_vs_virus_targets_db_all_hits.tsv \\
        \${nxf_workdir}/${meta.id}_all_hits.tsv
    cd \${nxf_workdir}

    # clean up intermediate mapping results to save disk
    rm -rf ${workflow.workDir.toAbsolutePath()}/spacerextractor_maptotarget/${meta.id}/
    """
}

process SPACEREXTRACTOR_COMBINEMAPS {
    label 'process_single'
    storeDir "."

    input:
    path(map_tsvs)

    output:
    path("${params.output}")    , emit: final_output

    script:
    """
    # iterate over phist tables
    for table in ${map_tsvs[0]}; do
       head -n 1 \${table} > ${params.output}
    done

    for table in ${map_tsvs}; do
        tail -n +2 \${table} >> ${params.output}
    done
    """
}


// Run entry workflow
workflow {
    main:
    // Check if output file already exists
    def output_file = file("${params.output}")

    if (!output_file.exists()) {

        // 1. Load/Download input files
        ch_virus_fasta = channel.fromPath(params.virus_fasta)
        ch_spacer_fasta = channel.fromPath(params.spacer_fasta)

        // 2. Split spacer file into chunks (process - SEQKIT_SPLIT2)
        SEQKIT_SPLIT2(
            ch_spacer_fasta
        )

        ch_split_fastas = SEQKIT_SPLIT2.out.split_fastas
            .map { file -> file }
            .flatten()
            .map { file ->
                [ [ id: file.getBaseName() ], file ]
            }

        // 3. Create spacerextractor DB from virus fasta (process - SPACEREXTRACTOR_CREATETARGETDB)
        SPACEREXTRACTOR_CREATETARGETDB(
            ch_virus_fasta
        )

        // 4. Align CRISPR chunks to spacerextractor DB (process - SPACEREXTRACTOR_MAPTOTARGET)
        SPACEREXTRACTOR_MAPTOTARGET(
            ch_split_fastas,
            SPACEREXTRACTOR_CREATETARGETDB.out.db.collect()
        )

        // 5. Combine all mapping results into single output file (process - SPACEREXTRACTOR_COMBINEMAPS)
        SPACEREXTRACTOR_COMBINEMAPS(
            SPACEREXTRACTOR_MAPTOTARGET.out.mapped_results.map { _meta, tsvs -> [ tsvs ] }.collect()
        )

    } else {
        println "Output file [${params.output}] already exists! Skipping nf-crisprhost."
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


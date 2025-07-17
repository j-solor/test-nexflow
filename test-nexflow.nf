#!/usr/bin/env nextflow

/*
 * Use echo to print 'Hello World!' to a file
 */
process BuildDb {

    publishDir 'database', mode: 'copy'

    input:
        val database_file

    output:
     
    path "my_local_db.*"  // Outputs all DB-related files (nsq, nin, nhr, etc.)

    script:
    """
    makeblastdb -in $database_file -dbtype nucl -out my_local_db
    """
}

process RunBLAST{
      publishDir 'results', mode: 'copy'

    input:
        path query_file
        path db_files // no need to call them latterl just make them available

    output:
        path "results.txt"

    script:
    """
    blastn -query ${query_file} -db my_local_db -out results.txt -outfmt 6
    """
}

process GetInfo{
    publishDir 'results', mode: 'copy'

    input:
        path blast_results

    output:
        path "obtined3UTRs.txt"

    script:
    """
    cat $blast_results > obtined3UTRs.txt   
    """

}

/*
 * Pipeline parameters
 */
params.database_file = 'input/toydb.fa'
params.query = 'input/misterious_sequences.fa'

workflow {

    database_file = Channel.fromPath(params.database_file)
    query = Channel.fromPath(params.query)

    // Connect processes
    db_files = BuildDb(database_file) // db_files will include .nin, .nsq, .nhr

    blast_results = RunBLAST(query, db_files)

    GetInfo(blast_results)
}

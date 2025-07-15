#!/usr/bin/env nextflow

/*
 * Use echo to print 'Hello World!' to a file
 */
process BuildDb {

    publishDir 'database', mode: 'copy'

    input:
        val database_file

    output:
     
     path "my_local_db*"
    
    script:
    """
    makeblastdb -in $database_file -dbtype nucl -out my_local_db
    """
}

process RunBLAST{
      publishDir 'results', mode: 'copy'

    input:
        path query_file

    output:
        path "results.txt"

    script:
    """
    blastn -query ${query_file} -db database/my_local_db -out results.txt -outfmt 6
    """
}

/*
 * Pipeline parameters
 */
params.database_file = 'input/toydb.fa'
params.query = 'input/misterious_sequences.fa'

workflow {

    // create a channel for inputs from a CSV file
    database_file = Channel.fromPath(params.database_file)
    query = Channel.fromPath(params.query)
 

    // Build de db
    BuildDb(database_file) // no output needed

    // Run Blast on fq files
    RunBLAST(query)

}

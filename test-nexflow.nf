#!/usr/bin/env nextflow

/*
 * Use echo to print 'Hello World!' to a file
 */
process BuildDb {

    publishDir 'results', mode: 'copy'

    input:
        val database_file

    script:
    """
    makeblastdb -in $database_file -dbtype nucl -out my_local_db
    """
}

process RunBLAST{
      publishDir 'results', mode: 'copy'

    input:
        path input_file

    output:
        path "UPPER-${input_file}-output.txt"

    script:
    """
    cat $input_file | tr '[a-z]' '[A-Z]' > UPPER-${input_file}-output.txt
    """
}

/*
 * Pipeline parameters
 */
params.database_file = 'input/toydb.fa'

workflow {

    // create a channel for inputs from a CSV file
    database_file = Channel.fromPath(params.database_file)
 

    // emit a greeting
    BuildDb(database_file)
}

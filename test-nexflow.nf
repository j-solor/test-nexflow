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
        val  blast_results

    output:
        path "obtined3UTRs.txt"

    script:
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(biomaRt)

    print("hola")
    #load Ensembl
    ensembl95 <- useEnsembl(
    biomart = "genes",
    dataset = "hsapiens_gene_ensembl"
    ) 

    annot_ensembl95 <- getBM(attributes = c(
    "ensembl_gene_id",
    "external_gene_name"),
    mart = ensembl95)

    translate <- deframe(annot_ensembl95[c("ensembl_gene_id", "external_gene_name")])

    fileConn<-file("obtined3UTRs.txt")
    writeLines(translate[c('$blast_results')], fileConn)
    close(fileConn)


    #! Now just need to turn this into a R script with my Biomart code
    #! Check speed, maybe better keep all outputs together and run BiomaRt once
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
    BuildDb(database_file) // db_files will include .nin, .nsq, .nhr
    RunBLAST(query, BuildDb.out)
    ensemblIDs = RunBLAST.out
                        .view {it -> "1st OUTPUT: $it" }
                        .splitCsv(header: false, sep:'\t')
                        .map{row -> "${row[1]}"}
                        .view {it -> "2nd OUTPUT: $it" }

              
    GetInfo(ensemblIDs)
    //! add a message with some info
}
    
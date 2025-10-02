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
        path "*_obtined3UTRs.tsv"

    script:
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(biomaRt)

    #load Ensembl
    ensembl95 <- useEnsembl(
    biomart = "genes",
    dataset = "hsapiens_gene_ensembl"
    ) 

    #load Ensembl
    ensembl95 <- useEnsembl( # listAttributes(ensembl95)
    biomart = "genes",
     dataset = "hsapiens_gene_ensembl") 

    chosen_atr <- c(
     "ensembl_gene_id",
    "ensembl_transcript_id",
    "external_gene_name",
     "chromosome_name",
     "3_utr_start",
     "3_utr_end",
     "3utr")

    query_ensemblIDs <- c(${blast_results.collect{"\"${it}\""}.join(", ")})

    annot_ensembl95 <- getBM(filters = c("ensembl_gene_id"),
                         values = query_ensemblIDs,
                         attributes = chosen_atr,
                         mart = ensembl95) %>%
    group_by(ensembl_gene_id)


    tibble_list <- annot_ensembl95 %>%
    group_split(ensembl_gene_id) %>% 
    set_names(group_keys(annot_ensembl95, ensembl_gene_id) %>% deframe())


    # Write results to file
    for(key in group_keys(annot_ensembl95, ensembl_gene_id) %>% deframe()){
    write_tsv(x = tibble_list[[key]], file = paste0(key,"_obtined3UTRs.tsv"))
    }

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
                        .collect()
                        .view()

              
    GetInfo(ensemblIDs)
    //! add a message with some info
}
    
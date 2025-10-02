
# Nexflow first workflow
This is a toy example workflow relying on the usage of Nextflow. For this, the task of choice is to get some misterious sequences and (1) find the gene they belong to (if they do to any) via blast and (2) Explore the 3'UTR of these given genes, retrieving their sequence and position.

Everything is self contained in `test-nextflow.nf`. It consists of 3 processes

1) process BuildDb: Constructs a toy database where blast will look for matches later. In this case it only contains the genomic sequence of the two proteins we are supposed to find and is `input/toydb.fa`

2) process RunBLAST: The constructed database is queried for the sequences of choice, in this case located in `input/misterious_sequences.fa`

3) process GetInfo: The results from blast are automatically parsed into R, in order to query `biomaRt` and obtain the 3'UTRs of every ensemblID that matched that target. BiomaRts is queried just once, which is specially important to not get your IP banned.

## Further developments
While being a test, many more things could be implemented. Examples of this would be: 
- modular input, where `input/toydb.fa` could be selected in the command line.
- In a real world scenario another variable should breach to the final table containing the name of the original sequence.
- Would be interesting to query BLAST directly.
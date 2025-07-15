#!/bin/bash

# Query the db to search for homology
blastn -query input/misterious_sequences.fa -db my_local_db -out results/results.txt -outfmt 6
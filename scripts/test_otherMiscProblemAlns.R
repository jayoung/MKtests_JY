## what happens with a ragged alignment? (seqs are not all the same length)
library(here)
library(Biostrings)
source("scripts/MKfunctions.R")


raggedAlnFile <- here("test_data/test_otherMiscProblemAlns/test_raggedAln.fa")

raggedAln <- readBStringSet(raggedAlnFile)
names(raggedAln) <- sapply(strsplit(names(raggedAln), " "), "[[", 1)


## get IDs of individuals in each population
raggedAlnSpecies <- gsub("\\d","",names(raggedAln))
raggedAlnPopulations <- split(names(raggedAln), raggedAlnSpecies)
 

#### single codon aln, without Ns
MKresults <- doMKtest(raggedAlnFile, 
                      pop1seqs=raggedAlnPopulations[["pongo"]],
                      pop2seqs=raggedAlnPopulations[["trachy"]])


## what happens with an aln where seqs are not all an even multiple of 3 in length?
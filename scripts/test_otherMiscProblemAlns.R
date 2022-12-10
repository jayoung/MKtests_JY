library(here)
library(Biostrings)
source("scripts/MKfunctions.R")

##### what happens with a ragged alignment? (seqs are not all the same length)
## I have now built in an error check for this

raggedAlnFile <- here("test_data/test_otherMiscProblemAlns/test_raggedAln.fa")

## get IDs of individuals in each population
raggedAln <- readBStringSet(raggedAlnFile)
names(raggedAln) <- sapply(strsplit(names(raggedAln), " "), "[[", 1)
raggedAlnSpecies <- gsub("\\d","",names(raggedAln))
raggedAlnPopulations <- split(names(raggedAln), raggedAlnSpecies)
 
# test
MKresults <- doMKtest(raggedAlnFile, 
                      pop1seqs=raggedAlnPopulations[["pongo"]],
                      pop2seqs=raggedAlnPopulations[["trachy"]])


####### what happens with an aln where alignment is not an even multiple of 3 in length?  I already check for that.
partialCodonAlnFile <- here("test_data/test_otherMiscProblemAlns/test_partialCodonAln.fa")

MKresults <- doMKtest(partialCodonAlnFile, 
                      pop1seqs=raggedAlnPopulations[["pongo"]],
                      pop2seqs=raggedAlnPopulations[["trachy"]])


####### mixed upper and lower case seqs
upperLowerAlnFile <- here("test_data/test_otherMiscProblemAlns/MKTwebsite_testAln_mixUpperLowerCase.fa")

# get IDs of individuals in each population
upperLowerAln <- readBStringSet(upperLowerAlnFile)
names(upperLowerAln) <- sapply(strsplit(names(upperLowerAln), " "), "[[", 1)
upperLowerAlnSpecies <- gsub("\\d","",names(upperLowerAln))
upperLowerAlnPopulations <- split(names(upperLowerAln), upperLowerAlnSpecies)

# test
MKresults <- doMKtest(upperLowerAlnFile, 
                      pop1seqs=upperLowerAlnPopulations[["pongo"]],
                      pop2seqs=upperLowerAlnPopulations[["trachy"]],
                      polarize = TRUE, outgroupSeqs = list(upperLowerAlnPopulations[["trachy"]]))


######### user might accidentally specify populations with names in common. that'll mess things up.
MKresults <- doMKtest(upperLowerAlnFile, 
                      pop1seqs=upperLowerAlnPopulations[["pongo"]],
                      pop2seqs=upperLowerAlnPopulations[["pongo"]])

####### unusual - an empty alignment file - I now handle this
emptyAlnFile <- here("test_data/test_otherMiscProblemAlns/test_emptyAlnFile.fa")

MKresults <- doMKtest(emptyAlnFile, 
                      pop1seqs=raggedAlnPopulations[["pongo"]],
                      pop2seqs=raggedAlnPopulations[["trachy"]])


######## stop codons fixed in one population but not the other (or in the outgroup)
stopCodonAllSeqs <- here("test_data/test_otherMiscProblemAlns/MKTwebsite_testAln.fa")
stopCodonAln1 <- here("test_data/test_otherMiscProblemAlns/test_onePopHasStop.fa")
stopCodonAln2 <- here("test_data/test_otherMiscProblemAlns/test_outgroupHasStop.fa")

# test
MKresults <- doMKtest(stopCodonAllSeqs, 
                      pop1seqs=raggedAlnPopulations[["pongo"]],
                      pop2seqs=raggedAlnPopulations[["trachy"]])


# test
MKresults <- doMKtest(stopCodonAln1, 
                      pop1seqs=raggedAlnPopulations[["pongo"]],
                      pop2seqs=raggedAlnPopulations[["trachy"]])




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



######## I want to be able to combine polarized and unpolarized results into the same output file


risaFiles <- list()
risaFiles[["popDat_justSim"]] <- here("data/Risa/2022_Dec7/abo_simpop_mel_yak.fasta")
risaFiles[["popDat_bothSpecies"]] <- here("data/Risa/2022_Dec7/abo_melpop_simpop_yak.fasta")

risa_alns <- lapply(risaFiles, readDNAStringSet)

## set up population names
risa_popIDs <- list()

risa_popIDs[["popDat_justSim"]] <- list()
risa_popIDs[["popDat_justSim"]][["sim"]] <- grep("CM015606", names(risa_alns[["popDat_justSim"]]), value=TRUE)
risa_popIDs[["popDat_justSim"]][["mel"]] <- grep("melanogaster", names(risa_alns[["popDat_justSim"]]), value=TRUE)
risa_popIDs[["popDat_justSim"]][["yak"]] <- grep("yakuba", names(risa_alns[["popDat_justSim"]]), value=TRUE)

risa_popIDs[["popDat_bothSpecies"]] <- list()
risa_popIDs[["popDat_bothSpecies"]][["sim"]] <- grep("CM015606|sim", names(risa_alns[["popDat_bothSpecies"]]), value=TRUE)
risa_popIDs[["popDat_bothSpecies"]][["yak"]] <- grep("yakuba", names(risa_alns[["popDat_bothSpecies"]]), value=TRUE)
## mel is everything else:
risa_popIDs[["popDat_bothSpecies"]][["mel"]] <- setdiff(names(risa_alns[["popDat_bothSpecies"]]), 
                                                        c(risa_popIDs[["popDat_bothSpecies"]][["sim"]],
                                                          risa_popIDs[["popDat_bothSpecies"]][["yak"]])) 
## sanity check
# lapply(risa_popIDs, function(x) { sapply(x, length)})

## MK tests
risa_MKresults <- list()

## popDat_justSim
risa_MKresults[["popDat_justSim_unpolarizedA"]] <- doMKtest(
    risaFiles[["popDat_justSim"]], 
    pop1seqs=risa_popIDs[["popDat_justSim"]][["sim"]], pop1alias="sim",
    pop2seqs=risa_popIDs[["popDat_justSim"]][["mel"]], pop2alias="mel")

risa_MKresults[["popDat_justSim_unpolarizedB"]] <- doMKtest(
    risaFiles[["popDat_justSim"]], 
    pop1seqs=risa_popIDs[["popDat_justSim"]][["mel"]], pop1alias="mel",
    pop2seqs=risa_popIDs[["popDat_justSim"]][["sim"]], pop2alias="sim")

risa_MKresults[["popDat_justSim_polarizedA"]] <- doMKtest(
    risaFiles[["popDat_justSim"]], 
    pop1seqs=risa_popIDs[["popDat_justSim"]][["sim"]], pop1alias="sim",
    pop2seqs=risa_popIDs[["popDat_justSim"]][["mel"]], pop2alias="mel",
    polarize=TRUE, 
    outgroupSeqs = list(risa_popIDs[["popDat_justSim"]][["yak"]]))

risa_MKresults[["popDat_justSim_polarizedB"]] <- doMKtest(
    risaFiles[["popDat_justSim"]], 
    pop1seqs=risa_popIDs[["popDat_justSim"]][["mel"]], pop1alias="mel",
    pop2seqs=risa_popIDs[["popDat_justSim"]][["sim"]], pop2alias="sim",
    polarize=TRUE, 
    outgroupSeqs = list(risa_popIDs[["popDat_justSim"]][["yak"]]))


## popDat_bothSpecies
risa_MKresults[["popDat_bothSpecies_unpolarizedA"]] <- doMKtest(
    risaFiles[["popDat_bothSpecies"]], 
    pop1seqs=risa_popIDs[["popDat_bothSpecies"]][["sim"]], pop1alias="sim",
    pop2seqs=risa_popIDs[["popDat_bothSpecies"]][["mel"]], pop2alias="mel")

risa_MKresults[["popDat_bothSpecies_unpolarizedB"]] <- doMKtest(
    risaFiles[["popDat_bothSpecies"]], 
    pop1seqs=risa_popIDs[["popDat_bothSpecies"]][["mel"]], pop1alias="mel",
    pop2seqs=risa_popIDs[["popDat_bothSpecies"]][["sim"]], pop2alias="sim")

risa_MKresults[["popDat_bothSpecies_polarizedA"]] <- doMKtest(
    risaFiles[["popDat_bothSpecies"]], 
    pop1seqs=risa_popIDs[["popDat_bothSpecies"]][["sim"]], pop1alias="sim",
    pop2seqs=risa_popIDs[["popDat_bothSpecies"]][["mel"]], pop2alias="mel",
    polarize=TRUE, 
    outgroupSeqs = list(risa_popIDs[["popDat_bothSpecies"]][["yak"]]))

risa_MKresults[["popDat_bothSpecies_polarizedB"]] <- doMKtest(
    risaFiles[["popDat_bothSpecies"]], 
    pop1seqs=risa_popIDs[["popDat_bothSpecies"]][["mel"]], pop1alias="mel",
    pop2seqs=risa_popIDs[["popDat_bothSpecies"]][["sim"]], pop2alias="sim",
    polarize=TRUE, 
    outgroupSeqs = list(risa_popIDs[["popDat_bothSpecies"]][["yak"]]))


## combine and save to Excel files

## without excel file
risa_MKresults_unpolarized <- combineMKresults(risa_MKresults[ grep("_unpolarized", names(risa_MKresults)) ])
risa_MKresults_polarized <- combineMKresults(risa_MKresults[ grep("_polarized", names(risa_MKresults)) ])


### with excel file
risa_MKresults_unpolarized <- combineMKresults(risa_MKresults[ grep("_unpolarized", names(risa_MKresults)) ],
                                               outFile="risa_MKresults_unpolarized_v1.xlsx",
                                               outDir=here("test_data/test_otherMiscProblemAlns"))
### with excel file keepNA=FALSE, NAcharacter="", 
risa_MKresults_unpolarized <- combineMKresults(risa_MKresults[ grep("_unpolarized", names(risa_MKresults)) ],
                                               outFile="risa_MKresults_unpolarized_v2.xlsx",
                                               outDir=here("test_data/test_otherMiscProblemAlns"),
                                               keepNA=FALSE, NAcharacter="")


risa_MKresults_all <- combineMKresults(risa_MKresults, extraVerbose=TRUE)
                                   
risa_MKresults_all <- combineMKresults(risa_MKresults, 
                                       outFile="risa_MKresults_all.xlsx",
                                       outDir=here("test_data/test_otherMiscProblemAlns"), 
                                       extraVerbose=TRUE)


### weird - ALL OF the column headers are blank with this:
risa_MKresults_all_2 <- combineMKresults(risa_MKresults, 
                                       outFile="risa_MKresults_all_v2.xlsx",
                                       outDir=here("test_data/test_otherMiscProblemAlns"),
                                       keepNA=FALSE, NAcharacter="", 
                                       extraVerbose=TRUE)

### weird - ALL OF the column headers are are blank using this:
risa_MKresults_all_3 <- combineMKresults(risa_MKresults, 
                                         outFile="risa_MKresults_all_v3.xlsx",
                                         outDir=here("test_data/test_otherMiscProblemAlns"),
                                         keepNA=FALSE, NAcharacter="N.A.", 
                                         extraVerbose=TRUE)





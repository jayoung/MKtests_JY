### https://en.wikipedia.org/wiki/McDonald%E2%80%93Kreitman_test
# alpha represents the proportion of substitutions driven by positive selection.
# neutrality index: >1 indicates negative selection; <1 indicates positive selection

### set working directory
# mac location
setwd("~/Desktop/mac_workStuff/mac_MKtests/MKtests_JY/")
# fast location
#setwd("/fh/fast/malik_h/user/jayoung/MKtest/forGithub/")

### load functions
source("scripts/MKfunctions.R")



##### unpolarized tests on Rini files
Rini_alnDir <- "data/Rini/alignments"
Rini_alnFiles <- list.files(Rini_alnDir, pattern="noStops.degapcodon.fas$", full.names=TRUE)
names(Rini_alnFiles) <- gsub(paste(Rini_alnDir, "/", sep=""), "", Rini_alnFiles)
names(Rini_alnFiles) <- gsub(".plusOutgroup_aln\\d_NT.edit.+?noStops.degapcodon.fas$", "", names(Rini_alnFiles))
names(Rini_alnFiles) <- gsub(".JYaln2.edit.noStops.degapcodon.fas$", "", names(Rini_alnFiles))

Rini_alns <- lapply(Rini_alnFiles, function(x) {
    aln <- readBStringSet(x)
    names(aln) <- sapply(strsplit(names(aln)," "),"[[",1)
    return(aln)
})

# get names of seqs within each population
Rini_populationNames <- lapply(Rini_alns, function(x){
    populationNames <- list()
    populationNames[["pop2"]] <- grep("sim|^GD", names(x), value=TRUE)
    outgroups <- grep("yak|ere", names(x), value=TRUE)
    populationNames[["out"]] <- list(outgroups)
    # all the rest are D. mel
    populationNames[["pop1"]] <- setdiff(names(x), 
                                         c(populationNames[["pop2"]], outgroups) )
    return(populationNames)
})
# sapply(Rini_populationNames, function(x) { sapply(x, length)})


## unpolarized test on all Rini alignments - works, gives same output as before
MKresults_RiniAlignments <- lapply (names(Rini_alnFiles), function(x){
    doMKtest(Rini_alnFiles[x], 
             outDir="data/Rini/results_unpolarized",
             pop1seqs=Rini_populationNames[[x]][["pop1"]], pop1alias="mel",
             pop2seqs=Rini_populationNames[[x]][["pop2"]], pop2alias="sim" )
})
names(MKresults_RiniAlignments) <- names(Rini_alnFiles)

## combine unpolarized results from all Rini alignments
Rini_results_all <- combineMKresults(MKresults_RiniAlignments, 
                                     outDir="data/Rini/results_unpolarized",
                                     outFile="Rini_allMKresults_unpolarized.xlsx", 
                                     getGeneNames=TRUE, 
                                     geneNameFile="data/Rini/riniTable2_geneOrder.txt",
                                     pop1alias="mel", pop2alias="sim")

## polarized test on all Rini alignments - works, gives same output as before
MKresults_RiniAlignments_polarized <- lapply (names(Rini_alnFiles), function(x){
    doMKtest(Rini_alnFiles[x], 
             outDir="data/Rini/results_polarized",
             pop1seqs=Rini_populationNames[[x]][["pop1"]], pop1alias="mel",
             pop2seqs=Rini_populationNames[[x]][["pop2"]], pop2alias="sim",
             outgroupSeqs=Rini_populationNames[[x]][["out"]],
             polarize=TRUE )
})
names(MKresults_RiniAlignments_polarized) <- names(Rini_alnFiles)

## combine polarized results from all Rini alignments
Rini_results_polarized_all <- combineMKresults(MKresults_RiniAlignments_polarized, 
                                               outDir="data/Rini/results_polarized",
                                               outFile="Rini_allMKresults_polarized.xlsx", 
                                               getGeneNames=TRUE, 
                                               geneNameFile="data/Rini/riniTable2_geneOrder.txt",
                                               pop1alias="mel", pop2alias="sim")

### https://en.wikipedia.org/wiki/McDonald%E2%80%93Kreitman_test
# alpha represents the proportion of substitutions driven by positive selection.
# neutrality index: >1 indicates negative selection; <1 indicates positive selection

### set working directory
# mac location - local
#setwd("~/Desktop/mac_workStuff/mac_MKtests/MKtests_JY")
# mac location - server
setwd("/Volumes/malik_h/user/jayoung/MKtest/MKtests_JY")
# fast location
#setwd("/fh/fast/malik_h/user/jayoung/MKtest/MKtests_JY")

### load functions
source("scripts/MKfunctions.R")


##### Gustavo's alignments, unpolarized tests

## v2 of alignments, see email Aug 3 2021.
# we have:
# 11 D. buzzatii populations (was 5)
# 3  D. seriema (as before changed) 
# outgroup species (D. stalkeri).


### now do MK tests
gustavo_alnDir <- "data/Gustavo/v2_2021_Aug16/originalFilesAsFasta"

gustavo_alnFiles <- list.files(gustavo_alnDir, pattern=".fa$", full.names=TRUE)
names(gustavo_alnFiles) <- gsub(paste(gustavo_alnDir, "/", sep=""), "", gustavo_alnFiles)
names(gustavo_alnFiles) <- gsub(".fa$", "", names(gustavo_alnFiles))
names(gustavo_alnFiles) <- gsub("nogap", "", names(gustavo_alnFiles))
names(gustavo_alnFiles) <- gsub("outgroup", "", names(gustavo_alnFiles))
names(gustavo_alnFiles) <- gsub("Cid6", "Cid6_", names(gustavo_alnFiles))
names(gustavo_alnFiles) <- gsub("Cid5", "Cid5_", names(gustavo_alnFiles))

gustavo_alns <- lapply(gustavo_alnFiles, function(x) {
    aln <- readBStringSet(x)
    names(aln) <- sapply(strsplit(names(aln)," "),"[[",1)
    return(aln)
})

# get names of seqs within each population
gustavo_populationNames <- lapply(gustavo_alns, function(x){
    populationNames <- list()
    populationNames[["pop1"]] <- grep("buz", names(x), value=TRUE)
    populationNames[["pop2"]] <- grep("sma", names(x), value=TRUE)
    populationNames[["out"]] <- list()
    populationNames[["out"]][[1]] <- grep("Dstl", names(x), value=TRUE)
    return(populationNames)
})
# sapply(gustavo_populationNames, function(x) { sapply(x, length)})

## unpolarized test on all Gustavo alignments - works
gustavo_results <- lapply (names(gustavo_alnFiles), function(x){
    doMKtest(gustavo_alnFiles[x], 
             outDir="data/Gustavo/v2_2021_Aug16/MKresults_unpolarized",
             pop1seqs=gustavo_populationNames[[x]][["pop1"]], pop1alias="buz",
             pop2seqs=gustavo_populationNames[[x]][["pop2"]], pop2alias="sma" )
})
names(gustavo_results) <- names(gustavo_alnFiles)

## combine results from all Gustavo alignments
gustavo_results_all <- combineMKresults(gustavo_results, 
                                        outFile="gustavo_allMKresults_unpolarized.xlsx",
                                        outDir="data/Gustavo/v2_2021_Aug16/MKresults_unpolarized",
                                        getGeneNames=FALSE, 
                                        pop1alias="buz", pop2alias="sma")

## polarized test on all Gustavo alignments - works
gustavo_results_polarized <- lapply (names(gustavo_alnFiles), function(x){
    doMKtest(gustavo_alnFiles[x], 
             outDir="data/Gustavo/v2_2021_Aug16/MKresults_polarized",
             pop1seqs=gustavo_populationNames[[x]][["pop1"]], pop1alias="buz",
             pop2seqs=gustavo_populationNames[[x]][["pop2"]], pop2alias="sma",
             outgroupSeqs=gustavo_populationNames[[x]][["out"]],
             polarize=TRUE )
})
names(gustavo_results_polarized) <- names(gustavo_alnFiles)

## combine results from all Gustavo alignments
gustavo_results_polarized_all <- combineMKresults(gustavo_results_polarized, 
                                                  outFile="gustavo_allMKresults_polarized.xlsx",
                                                  outDir="data/Gustavo/v2_2021_Aug16/MKresults_polarized",
                                                  getGeneNames=FALSE, 
                                                  pop1alias="buz", pop2alias="sma")


## xx not plotting!
pdf(height=5,width=15,file="data/Gustavo/v2_2021_Aug16/MKresults_polarized/gustavo_allMKresults_polarized.pdf")
## using layout rather than par-mfrow so that plots get added in columns before rows
layout(matrix(1:9, nrow=3, ncol=3))
for (x in names(gustavo_results_polarized)) {
    cat("plotting",x,"\n")
    plotMKpositions(gustavo_results_polarized[[x]][["positions"]], 
                    plotPolarizedPop1=TRUE, plotPolarizedPop2=TRUE, 
                    title=x, pop1alias="buz", pop2alias="sma", 
                    setNumPlots=FALSE)
}
dev.off()



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


##### chingHo's alignment, unpolarized tests

chingHo_alnFile <- "data/ChingHo/BigH1_CDS_with_sim.fasta"
chingHo_aln <- readBStringSet(chingHo_alnFile)

numSeqs <- length(chingHo_aln)

## get names of seqs for species 1 (pop1 = D. melanogaster) and species 2 (pop2 = D. simulans, the last sequence)
chingHo_populationNames <- list()
chingHo_populationNames[["pop1"]] <- names(chingHo_aln)[1:(numSeqs-1)]
chingHo_populationNames[["pop2"]] <- names(chingHo_aln)[numSeqs]

## simple unpolarized test on all chingHo alignments 
chingHo_results <- doMKtest(chingHo_alnFile, 
             outDir="data/ChingHo/MKresults_unpolarized",
             pop1seqs=chingHo_populationNames[["pop1"]], pop1alias="Dmel",
             pop2seqs=chingHo_populationNames[["pop2"]], pop2alias="Dsim" )


## same, but remove rare polymorphisms (<=5%)
chingHo_results_noRare <- doMKtest(chingHo_alnFile, 
                            outDir="data/ChingHo/MKresults_unpolarized",
                            pop1seqs=chingHo_populationNames[["pop1"]], pop1alias="Dmel",
                            pop2seqs=chingHo_populationNames[["pop2"]], pop2alias="Dsim", 
                            filterRareAlleles=TRUE, alleleFreqThreshold=0.05)



####### perhaps we want to do it on single D. mel populations at once, rather than individuals from a mix of populations:

## first we get population of each strain from its name (exceptions:  BigH1 is the reference Dmel sequence). This will be a character vector in the same order as chingHo_populationNames[["pop1"]], so we can use it to select sequences we care about
populationNames <- sapply(strsplit(chingHo_populationNames[["pop1"]], "_"), "[[", 1) 
populationNames <- gsub("\\d.+?$","",populationNames)
populationNames <- gsub("-$","",populationNames)

## how many individuals in each population?
table(populationNames)

## get populations with at least 20 individuals:
bigEnoughPopulations <- names(which( table(populationNames) >= 20 ))

## now get the names of strains for each of those populations:
bigEnoughPopulations_strainNames <- lapply( bigEnoughPopulations, function(x) {
    chingHo_populationNames[["pop1"]][which(populationNames==x)]
})
names(bigEnoughPopulations_strainNames) <- bigEnoughPopulations



## unpolarized test for each of those populations
chingHo_resultsEachPopulation <- lapply(names(bigEnoughPopulations_strainNames), function(popName) {
    cat("doing MK test for",popName,"population\n")
    doMKtest(chingHo_alnFile, 
             outDir="data/ChingHo/MKresults_unpolarized",
             pop1seqs=bigEnoughPopulations_strainNames[[popName]], pop1alias=paste("Dmel",popName,sep="_"),
             pop2seqs=chingHo_populationNames[["pop2"]], pop2alias="Dsim" )
})
names(chingHo_resultsEachPopulation) <- paste("BigH1",names(bigEnoughPopulations_strainNames),sep="_")


## unpolarized test for each of those populations removing rare variants (<=5%)
chingHo_resultsEachPopulation_noRare <- lapply(names(bigEnoughPopulations_strainNames), function(popName) {
    cat("doing MK test for",popName,"population\n")
    doMKtest(chingHo_alnFile, 
             outDir="data/ChingHo/MKresults_unpolarized",
             pop1seqs=bigEnoughPopulations_strainNames[[popName]], pop1alias=paste("Dmel",popName,sep="_"),
             pop2seqs=chingHo_populationNames[["pop2"]], pop2alias="Dsim", 
             filterRareAlleles=TRUE, alleleFreqThreshold=0.05 )
})
names(chingHo_resultsEachPopulation_noRare) <- paste("BigH1",names(bigEnoughPopulations_strainNames),"noRare",sep="_")


##### combine all those results together and save as an Excel file
allResults <- list(BigH1_all=chingHo_results, BigH1_all_noRare=chingHo_results_noRare)
allResults <- c(allResults, chingHo_resultsEachPopulation, chingHo_resultsEachPopulation_noRare)
combinedResults <- combineMKresults(allResults, 
                                    outFile="BigH1_MKtest_unpolarized_results.xlsx",
                                    outDir="data/ChingHo/MKresults_unpolarized",
                                    getGeneNames=FALSE, 
                                    pop1alias="Dmel", pop2alias="Dsim")


###### make plots that show where all the changes are for each result
pdf(height=7,width=11,file="data/ChingHo/MKresults_unpolarized/BigH1_MKtest_unpolarized_plots.pdf")        
par(mfrow=c(2,1))
for (x in c("BigH1_all_noRare", "BigH1_ZI_noRare")) {
    plotMKpositions(allResults[[x]][["positions"]],  
                    title=x, pop1alias="Dmel", pop2alias="Dsim", 
                    setNumPlots=FALSE)
}
dev.off()


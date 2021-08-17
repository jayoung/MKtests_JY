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


# Error in getCodonChangeCounts(codonsBefore, codonsAfter) : 
#     ERROR - found a codon missing from the paths table, probably contains a stop codon:GAGTAG at position 1  

# troubleshoot
# seqs 1-27 work OK, 
# seqs 1-28 does not work - give the codon missing from paths error

# bases 853-855 contain the problem
# seq 28 has a GAG codon

problemAln <- narrow(chingHo_aln, start=853, width=3)[c(1,16,28,1115)]
temp_results <- doMKtest(myAln=problemAln, 
                         outfileStem="tempMK", 
                         pop1seqs=chingHo_populationNames[["pop1"]][c(1,16,28)], pop1alias="Dmel",
                         pop2seqs=chingHo_populationNames[["pop2"]], pop2alias="Dsim", 
                         quiet=TRUE )

temp_results[[2]][,c("pop1_anc","pop2_anc","pop1_poly" ,"pop2_poly",
                     "pop1_vs_pop2_Dn", "pop1_vs_pop2_Ds", "pop1_Pn", "pop1_Ps")]



alnSlicesUniqueSeqs(problemAln)

writeXStringSet(problemAln, "problemAln.fa", format="fasta")

testAln <- BStringSet(c("GAG","TCG","GCG","GAG"))
names(testAln) <- c("seq1","seq2","seq3","seq4")
doMKtest(myAln=testAln, 
         outfileStem="tempMK", 
         pop1seqs=c("seq1","seq2"),
         pop2seqs=c("seq3"), pop2alias="Dsim", 
         quiet=TRUE )

## same, but remove rare polymorphisms (<=5%)
chingHo_results_noRare <- doMKtest(chingHo_alnFile, 
                            outDir="data/ChingHo/MKresults_unpolarized",
                            pop1seqs=chingHo_populationNames[["pop1"]], pop1alias="Dmel",
                            pop2seqs=chingHo_populationNames[["pop2"]], pop2alias="Dsim", 
                            filterRareAlleles=TRUE, alleleFreqThreshold=0.05)



####### perhaps we want to do it on single D. mel populations at once, rather than individuals from a mix of populations:

## first we get population of each strain from its name (exceptions:  Dp1 is the reference Dmel sequence). This will be a character vector in the same order as chingHo_populationNames[["pop1"]], so we can use it to select sequences we care about
populationNames <- sapply(strsplit(chingHo_populationNames[["pop1"]], "_"), "[[", 1) 
populationNames <- gsub("\\d.+?$","",populationNames)
populationNames <- gsub("-$","",populationNames)

## how many individuals in each population?
table(populationNames)
# CO CO1 CO2 Dp1  EA EA3  EB  ED ED2 ED3  EF EF2  EG  ER  EZ EZ2  FR  GA  GU GU2 GU6 GU7 GU9  KN KN6  KR KR7  NG NG7 NG9 RAL  RG RG2 
# 16   2   2   1  17   2   4   5   1   1  35   2  16   5   3   1  47  16   1   1   1   1   1   4   1   2   1   4   1   1  82  31   1 
# RG3 RG5 RG7 RG8 RG9  SB  SD  SF  SP  UG UG7  UK UK2  ZI  ZS ZS5 
#   1   1   1   1   1   5  24   3  28   3   1   3   1 180   3   1 

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
names(chingHo_resultsEachPopulation) <- paste("Dp1",names(bigEnoughPopulations_strainNames),sep="_")


## unpolarized test for each of those populations removing rare variants (<=5%)
chingHo_resultsEachPopulation_noRare <- lapply(names(bigEnoughPopulations_strainNames), function(popName) {
    cat("doing MK test for",popName,"population\n")
    doMKtest(chingHo_alnFile, 
             outDir="data/ChingHo/MKresults_unpolarized",
             pop1seqs=bigEnoughPopulations_strainNames[[popName]], pop1alias=paste("Dmel",popName,sep="_"),
             pop2seqs=chingHo_populationNames[["pop2"]], pop2alias="Dsim", 
             filterRareAlleles=TRUE, alleleFreqThreshold=0.05 )
})
names(chingHo_resultsEachPopulation_noRare) <- paste("Dp1",names(bigEnoughPopulations_strainNames),"noRare",sep="_")


##### combine all those results together and save as an Excel file
allResults <- list(Dp1_all=chingHo_results, Dp1_all_noRare=chingHo_results_noRare)
allResults <- c(allResults, chingHo_resultsEachPopulation, chingHo_resultsEachPopulation_noRare)
combinedResults <- combineMKresults(allResults, 
                                    outFile="Dp1_MKtest_unpolarized_results.xlsx",
                                    outDir="data/ChingHo/MKresults_unpolarized",
                                    getGeneNames=FALSE, 
                                    pop1alias="Dmel", pop2alias="Dsim")


###### make plots that show where all the changes are for each result
pdf(height=15,width=11,file="data/ChingHo/MKresults_unpolarized/Dp1_MKtest_unpolarized_plots.pdf")
for (x in names(allResults)) {
    cat("plotting",x,"\n")
    if (x %in% c("Dp1_all", "Dp1_EF", "Dp1_EF_noRare")) { 
        cat("    starting new page\n")
        par(mfrow=c(8,1))
    }
    plotMKpositions(allResults[[x]][["positions"]],  
                    title=x, pop1alias="Dmel", pop2alias="Dsim", 
                    setNumPlots=FALSE)
}
dev.off()


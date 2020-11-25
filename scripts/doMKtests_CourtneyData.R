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


######## run on three Arp53D alignments
Arp53D_alnDir <- "data/Arp53D"
Arp53D_alnFiles <- c(
    "CG5409-PB_Arp53D_popfly_strains_dsim_dyak_alignment_for_MKT.fasta.noStops.degapcodon",  
    "CG5409-PB_Arp53D_popfly_strains_dsim_dyak_alignment_for_MKT.fasta.noStops.degapcodon.trim118_135.fa", 
    "CG5409-PB_Arp53D_popfly_strains_dsim_dyak_alignment_for_MKT.fasta.noStops.degapcodon.justMel.noN.addSimYak.fa")
Arp53D_alnFiles <- paste(Arp53D_alnDir, Arp53D_alnFiles, sep="/")
names(Arp53D_alnFiles) <-c("fullLen", "justCodons40to45", "fullLenJust62mel")

Arp53D_alns <- lapply(Arp53D_alnFiles, function(x) {
    aln <- readBStringSet(x)
    names(aln) <- sapply(strsplit(names(aln)," "),"[[",1)
    return(aln)
})

# get names of seqs within each population
Arp53D_populationNames <- lapply(Arp53D_alns, function(x){
    populationNames <- list()
    populationNames[["pop1"]] <- grep("ZI|mel", names(x), value=TRUE)
    populationNames[["pop2"]] <- grep("sim", names(x), value=TRUE)
    populationNames[["out"]] <- grep("yak", names(x), value=TRUE)
    return(populationNames)
})
# sapply(Arp53D_populationNames, function(x) { sapply(x, length)})

### add another set to analyze, which is the noN version but without mel
Arp53D_alnFiles[["noN_justZI"]] <- Arp53D_alnFiles[["fullLenJust62mel"]]
Arp53D_alns[["noN_justZI"]] <- Arp53D_alns[["fullLenJust62mel"]]
Arp53D_populationNames[["noN_justZI"]] <- list()
Arp53D_populationNames[["noN_justZI"]][["pop1"]] <- grep("ZI", names(Arp53D_alns[["fullLenJust62mel"]]), value=TRUE)
Arp53D_populationNames[["noN_justZI"]][["pop2"]] <- grep("sim", names(Arp53D_alns[["fullLenJust62mel"]]), value=TRUE)
Arp53D_populationNames[["noN_justZI"]][["out"]] <- grep("yak", names(Arp53D_alns[["fullLenJust62mel"]]), value=TRUE)


### add another set to analyze, which is the full length version but without mel
Arp53D_alnFiles[["inclN_justZI"]] <- Arp53D_alnFiles[["fullLen"]]
Arp53D_alns[["inclN_justZI"]] <- Arp53D_alns[["fullLen"]]
Arp53D_populationNames[["inclN_justZI"]] <- list()
Arp53D_populationNames[["inclN_justZI"]][["pop1"]] <- grep("ZI", names(Arp53D_alns[["inclN_justZI"]]), value=TRUE)
Arp53D_populationNames[["inclN_justZI"]][["pop2"]] <- grep("sim", names(Arp53D_alns[["inclN_justZI"]]), value=TRUE)
Arp53D_populationNames[["inclN_justZI"]][["out"]] <- grep("yak", names(Arp53D_alns[["inclN_justZI"]]), value=TRUE)

### do various flavors of MK test on those three Arp53D alignments. I won't run the combineMK script on them, but instead I'll add the results to those of all Arps and get combined output

# unpolarized on all Arp53D datasets - works
Arp53D_allResults <- lapply (names(Arp53D_alnFiles), function(x){
    doMKtest(Arp53D_alnFiles[x], 
             pop1seqs=Arp53D_populationNames[[x]][["pop1"]], pop1alias="mel",
             pop2seqs=Arp53D_populationNames[[x]][["pop2"]], pop2alias="sim" )
})
names(Arp53D_allResults) <- names(Arp53D_alnFiles)

# unpolarized on all Arp53D datasets, with 5% freq filter
Arp53D_freqFilt_5pc_allResults <- lapply (names(Arp53D_alnFiles), function(x){
    doMKtest(Arp53D_alnFiles[x], 
             pop1seqs=Arp53D_populationNames[[x]][["pop1"]], pop1alias="mel",
             pop2seqs=Arp53D_populationNames[[x]][["pop2"]], pop2alias="sim",
             filterRareAlleles=TRUE, alleleFreqThreshold=0.05 )
})
names(Arp53D_freqFilt_5pc_allResults) <- names(Arp53D_alnFiles)

# polarized on all Arp53D datasets, with 5% freq filter
Arp53D_polarize_freqFilt_5pc_allResults <- lapply (names(Arp53D_alnFiles), function(x){
    doMKtest(Arp53D_alnFiles[x], 
             pop1seqs=Arp53D_populationNames[[x]][["pop1"]], pop1alias="mel",
             pop2seqs=Arp53D_populationNames[[x]][["pop2"]], pop2alias="sim",
             polarize=TRUE, outgroupSeqs = list(Arp53D_populationNames[[x]][["out"]]),
             filterRareAlleles=TRUE, alleleFreqThreshold=0.05 )
})
names(Arp53D_polarize_freqFilt_5pc_allResults) <- names(Arp53D_alnFiles)



##### all Courtney Arp alignments, before removing seqs with Ns, unpolarized tests
Arps_inclNseqs_alnDir <- "data/Arps/afterFixingArp2/aligned/byPopulation/ZI"

Arps_inclNseqs_alnFiles <- list.files(Arps_inclNseqs_alnDir, pattern="ZI.fa$", full.names=TRUE)
names(Arps_inclNseqs_alnFiles) <- gsub(paste(Arps_inclNseqs_alnDir, "/", sep=""), "", Arps_inclNseqs_alnFiles)
names(Arps_inclNseqs_alnFiles) <- gsub("_minus_introns.plusSim.ZI.fa$", "", names(Arps_inclNseqs_alnFiles))

Arps_inclNseqs_alns <- lapply(Arps_inclNseqs_alnFiles, function(x) {
    aln <- readBStringSet(x)
    names(aln) <- sapply(strsplit(names(aln)," "),"[[",1)
    return(aln)
})

# get names of seqs within each population
Arps_inclNseqs_populationNames <- lapply(Arps_inclNseqs_alns, function(x){
    populationNames <- list()
    populationNames[["pop2"]] <- grep("Dsim", names(x), value=TRUE)
    # all the rest are D. mel
    populationNames[["pop1"]] <- setdiff(names(x), c(populationNames[["pop2"]], populationNames[["out"]]) )
    return(populationNames)
})
# sapply(Arps_inclNseqs_populationNames, function(x) { sapply(x, length)})

### various flavors of test on those alignments. Again, I won't run the combine scripts on these - instead I add the corresponding Arp53D result and combine

## unpolarized test, no freq filter
Arps_inclNseqs_results <- lapply (names(Arps_inclNseqs_alnFiles), function(x){
    doMKtest(Arps_inclNseqs_alnFiles[x], 
             pop1seqs=Arps_inclNseqs_populationNames[[x]][["pop1"]], pop1alias="mel",
             pop2seqs=Arps_inclNseqs_populationNames[[x]][["pop2"]], pop2alias="sim" )
})
names(Arps_inclNseqs_results) <- names(Arps_inclNseqs_alnFiles)

## unpolarized test, with 5% freq filter
Arps_inclNseqs_freqFilt_5pc_results <- lapply (names(Arps_inclNseqs_alnFiles), function(x){
    doMKtest(Arps_inclNseqs_alnFiles[x], 
             pop1seqs=Arps_inclNseqs_populationNames[[x]][["pop1"]], pop1alias="mel",
             pop2seqs=Arps_inclNseqs_populationNames[[x]][["pop2"]], pop2alias="sim",
             filterRareAlleles=TRUE, alleleFreqThreshold=0.05)
})
names(Arps_inclNseqs_freqFilt_5pc_results) <- names(Arps_inclNseqs_alnFiles)



##### all Courtney new alignments, after removing seqs with Ns, unpolarized tests
Arps_noNseqs_alnDir <- "data/Arps/afterFixingArp2/aligned/byPopulation/ZI/removeNseqs"

Arps_noNseqs_alnFiles <- list.files(Arps_noNseqs_alnDir, pattern="ZI.noNseqs.fa$", full.names=TRUE)
names(Arps_noNseqs_alnFiles) <- gsub(paste(Arps_noNseqs_alnDir, "/", sep=""), "", Arps_noNseqs_alnFiles)
names(Arps_noNseqs_alnFiles) <- gsub("_minus_introns.plusSim.ZI.noNseqs.fa$", "", names(Arps_noNseqs_alnFiles))

Arps_noNseqs_alns <- lapply(Arps_noNseqs_alnFiles, function(x) {
    aln <- readBStringSet(x)
    names(aln) <- sapply(strsplit(names(aln)," "),"[[",1)
    return(aln)
})

# get names of seqs within each population
Arps_noNseqs_populationNames <- lapply(Arps_noNseqs_alns, function(x){
    populationNames <- list()
    populationNames[["pop2"]] <- grep("Dsim", names(x), value=TRUE)
    # all the rest are D. mel
    populationNames[["pop1"]] <- setdiff(names(x), c(populationNames[["pop2"]], populationNames[["out"]]) )
    return(populationNames)
})
# sapply(Arps_noNseqs_populationNames, function(x) { sapply(x, length)})

### various flavors of tests. I won't combine the Arps, but will add the relevant Arp53D result and then combine
## unpolarized test on all Arps alignments - works
Arps_noNseqs_results <- lapply (names(Arps_noNseqs_alnFiles), function(x){
    doMKtest(Arps_noNseqs_alnFiles[x], 
             pop1seqs=Arps_noNseqs_populationNames[[x]][["pop1"]], pop1alias="mel",
             pop2seqs=Arps_noNseqs_populationNames[[x]][["pop2"]], pop2alias="sim" )
})
names(Arps_noNseqs_results) <- names(Arps_noNseqs_alnFiles)

## unpolarized test on all Arps alignments without N seqs, with 5% freq filter
Arps_noNseqs_freqFilt_5pc_results <- lapply (names(Arps_noNseqs_alnFiles), function(x){
    doMKtest(Arps_noNseqs_alnFiles[x], 
             pop1seqs=Arps_noNseqs_populationNames[[x]][["pop1"]], pop1alias="mel",
             pop2seqs=Arps_noNseqs_populationNames[[x]][["pop2"]], pop2alias="sim",
             filterRareAlleles=TRUE, alleleFreqThreshold=0.05)
})
names(Arps_noNseqs_freqFilt_5pc_results) <- names(Arps_noNseqs_alnFiles)



##### add the corresponding Arp53D results to those for the other Arps

# reorder:
Arps_preferredOrder <- c("Arp1", "Arp2", "Arp3", "Arp4", "Arp5", 
                         "Arp6", "Arp8", "Arp10", "Arp53D")

## noN seqs
ArpsPlusArp53D_noNseqs_results <- c( Arps_noNseqs_results, 
                                     Arp53D_allResults["noN_justZI"])
names(ArpsPlusArp53D_noNseqs_results) <- gsub("noN_justZI","Arp53D",names(ArpsPlusArp53D_noNseqs_results))
ArpsPlusArp53D_noNseqs_results <- ArpsPlusArp53D_noNseqs_results[Arps_preferredOrder]

ArpsPlusArp53D_noNseqs_results_all <- combineMKresults(ArpsPlusArp53D_noNseqs_results, 
                                                       outFile="Arps_plus_Arp53D_allMKresults_noNseqs_noFreqFilter.xlsx",
                                                       outDir="data/Arp53D",
                                                       getGeneNames=FALSE, 
                                                       pop1alias="mel", pop2alias="sim")

## inclN seqs
ArpsPlusArp53D_inclNseqs_results <- c( Arps_inclNseqs_results, Arp53D_allResults["inclN_justZI"])
names(ArpsPlusArp53D_inclNseqs_results) <- gsub("inclN_justZI","Arp53D",names(ArpsPlusArp53D_inclNseqs_results))
ArpsPlusArp53D_inclNseqs_results <- ArpsPlusArp53D_inclNseqs_results[Arps_preferredOrder]

ArpsPlusArp53D_inclNseqs_results_all <- combineMKresults(ArpsPlusArp53D_inclNseqs_results, 
                                                         outFile="Arps_plus_Arp53D_allMKresults_inclNseqs_noFreqFilter.xlsx",
                                                         outDir="data/Arp53D",
                                                         getGeneNames=FALSE, 
                                                         pop1alias="mel", pop2alias="sim")

## noN seqs, 5% freq filter
ArpsPlusArp53D_noNseqs_freqFilt_5pc_results <- c( Arps_noNseqs_freqFilt_5pc_results, 
                                                  Arp53D_freqFilt_5pc_allResults["noN_justZI"])
names(ArpsPlusArp53D_noNseqs_freqFilt_5pc_results) <- gsub("noN_justZI","Arp53D",names(ArpsPlusArp53D_noNseqs_freqFilt_5pc_results))
ArpsPlusArp53D_noNseqs_freqFilt_5pc_results <- ArpsPlusArp53D_noNseqs_freqFilt_5pc_results[Arps_preferredOrder]

ArpsPlusArp53D_noNseqs_freqFilt_5pc_results_all <- combineMKresults(
    ArpsPlusArp53D_noNseqs_freqFilt_5pc_results, 
    outFile="Arps_plus_Arp53D_allMKresults_noNseqs_freqFilter_5pc.xlsx",
    outDir="data/Arp53D",
    getGeneNames=FALSE, 
    pop1alias="mel", pop2alias="sim")

## incl N seqs, 5% freq filter
ArpsPlusArp53D_inclNseqs_freqFilt_5pc_results <- c( Arps_inclNseqs_freqFilt_5pc_results, 
                                                    Arp53D_freqFilt_5pc_allResults["inclN_justZI"])
names(ArpsPlusArp53D_inclNseqs_freqFilt_5pc_results) <- gsub("inclN_justZI","Arp53D",names(ArpsPlusArp53D_inclNseqs_freqFilt_5pc_results))

ArpsPlusArp53D_inclNseqs_freqFilt_5pc_results <- 
    ArpsPlusArp53D_inclNseqs_freqFilt_5pc_results[Arps_preferredOrder]

ArpsPlusArp53D_inclNseqs_freqFilt_5pc_results_all <- combineMKresults(
    ArpsPlusArp53D_inclNseqs_freqFilt_5pc_results, 
    outFile="Arps_plus_Arp53D_allMKresults_inclNseqs_freqFilter_5pc.xlsx",
    outDir="data/Arp53D",
    getGeneNames=FALSE, 
    pop1alias="mel", pop2alias="sim")



#### plots

### without freq filter
pdf(height=11,width=7, 
    file="data/Arp53D/ArpsPlusArp53D_inclNseqs_noFreqFilt_MKplots.pdf")
par(mfrow=c(5,1))
for (x in names(ArpsPlusArp53D_inclNseqs_results)) {
    cat("plotting",x,"\n")
    plotMKpositions(ArpsPlusArp53D_inclNseqs_results[[x]][["positions"]], 
                    title=x, 
                    setNumPlots=FALSE)
}
dev.off()

pdf(height=11,width=7, 
    file="data/Arp53D/ArpsPlusArp53D_noNseqs_noFreqFilt_MKplots.pdf")
par(mfrow=c(5,1))
for (x in names(ArpsPlusArp53D_noNseqs_results)) {
    cat("plotting",x,"\n")
    plotMKpositions(ArpsPlusArp53D_noNseqs_results[[x]][["positions"]], 
                    title=x, 
                    setNumPlots=FALSE)
}
dev.off()

### with freq filter
pdf(height=11,width=7, 
    file="data/Arp53D/ArpsPlusArp53D_inclNseqs_freqFilt_5pc_MKplots.pdf")
par(mfrow=c(5,1))
for (x in names(ArpsPlusArp53D_inclNseqs_freqFilt_5pc_results)) {
    cat("plotting",x,"\n")
    plotMKpositions(ArpsPlusArp53D_inclNseqs_freqFilt_5pc_results[[x]][["positions"]], 
                    title=x, 
                    setNumPlots=FALSE)
}
dev.off()

pdf(height=11,width=7, 
    file="data/Arp53D/ArpsPlusArp53D_noNseqs_freqFilt_5pc_MKplots.pdf")
par(mfrow=c(5,1))
for (x in names(ArpsPlusArp53D_noNseqs_freqFilt_5pc_results)) {
    cat("plotting",x,"\n")
    plotMKpositions(ArpsPlusArp53D_noNseqs_freqFilt_5pc_results[[x]][["positions"]], 
                    title=x, 
                    setNumPlots=FALSE)
}
dev.off()



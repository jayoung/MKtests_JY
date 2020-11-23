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

######## test using same alignment as MK website:
MKresults_websiteExample <- doMKtest("test_data/MKTwebsiteExample/MKTwebsite_testAln.fa", 
          pop1seqs=c("pongo1","pongo2"), pop1alias="pongo",
          pop2seqs=c("trachy1","trachy2"), pop2alias="trachy",
          writeAncFasta=TRUE, writeMKoutput=TRUE)
#temp[["summary"]]
#head(temp[["positions"]])

######## test some alignments I made to test rare variant filtering:
test1_file <- "test_data/test_rareVariants/test1/test1.fa"
test1_aln <- readBStringSet(test1_file)
test1_popNames <- list()
test1_popNames[["pop1"]] <- grep ("pop1", names(test1_aln), value=TRUE)
test1_popNames[["pop2"]] <- grep ("pop2", names(test1_aln), value=TRUE)

test1_results <- doMKtest(test1_file, pop1seqs=test1_popNames[["pop1"]],
                          pop2seqs=test1_popNames[["pop2"]],
                          writeAncFasta=TRUE, writeMKoutput=TRUE,
                          flagRareAlleles=TRUE, alleleFreqThreshold=0.2)


test1_results_freqFilter_20pc <- doMKtest(test1_file, pop1seqs=test1_popNames[["pop1"]],
                          pop2seqs=test1_popNames[["pop2"]],
                          writeAncFasta=TRUE, writeMKoutput=TRUE,
                          flagRareAlleles=TRUE, filterRareAlleles=TRUE,
                          alleleFreqThreshold=0.2)

#test1_results[["summary"]]
#head(test1_results[["positions"]])

### test2 
test2_file <- "test_data/test_rareVariants/test2/test2.fa"
test2_aln <- readBStringSet(test2_file)
test2_popNames <- list()
test2_popNames[["pop1"]] <- grep ("pop1", names(test2_aln), value=TRUE)
test2_popNames[["pop2"]] <- grep ("pop2", names(test2_aln), value=TRUE)

test2_results_noFreqFilter <- doMKtest(test2_file, pop1seqs=test2_popNames[["pop1"]],
                          pop2seqs=test2_popNames[["pop2"]],
                          writeAncFasta=TRUE, writeMKoutput=TRUE,
                          flagRareAlleles=TRUE, 
                          alleleFreqThreshold=0.2)

test2_results_freqFilter_20pc <- doMKtest(test2_file, pop1seqs=test2_popNames[["pop1"]],
                          pop2seqs=test2_popNames[["pop2"]],
                          writeAncFasta=TRUE, writeMKoutput=TRUE,
                          flagRareAlleles=TRUE, filterRareAlleles=TRUE,
                          alleleFreqThreshold=0.2)


### test3 (three alignments)
test3_dir <- "test_data/test_rareVariants/test3"
test3_files <- list.files(test3_dir, pattern=".fa$", full.names=TRUE)
names(test3_files) <- gsub(paste(test3_dir,"/",sep=""),"",test3_files)
names(test3_files) <- gsub(".fa$","",names(test3_files))

test3_alns <- lapply(test3_files, readBStringSet)

test3_popNames <- lapply(test3_alns, function(x){
    popNames <- list()
    popNames[["pop1"]] <- grep ("pop1", names(x), value=TRUE)
    popNames[["pop2"]] <- grep ("pop2", names(x), value=TRUE)
    return(popNames)
})

test3_results <- lapply(names(test3_files), function(x){
    doMKtest(test3_files[[x]], pop1seqs=test3_popNames[[x]][["pop1"]],
             pop2seqs=test3_popNames[[x]][["pop2"]],
             writeAncFasta=TRUE, writeMKoutput=TRUE,
             flagRareAlleles=TRUE, filterRareAlleles=TRUE, 
             alleleFreqThreshold=0.2)
})



## combine results from all Arps alignments
test3_results_all <- combineMKresults(test3_results, 
                                       outFile="test3_results_freqFilter_20pc.xlsx",
                                       outDir=test3_dir,
                                       getGeneNames=FALSE)



#test3_results[["summary"]]
#head(test3_results[["positions"]])


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


# unpolarized on all Arp53D datasets - works
Arp53D_allResults <- lapply (names(Arp53D_alnFiles), function(x){
    doMKtest(Arp53D_alnFiles[x], 
          pop1seqs=Arp53D_populationNames[[x]][["pop1"]], pop1alias="mel",
          pop2seqs=Arp53D_populationNames[[x]][["pop2"]], pop2alias="sim" )
})
names(Arp53D_allResults) <- names(Arp53D_alnFiles)

## combine results from all Arps alignments
Arp53D_results_all <- combineMKresults(Arp53D_allResults, 
                outFile="Arp53D_allMKresults_conservative.xlsx",
                outDir=Arp53D_alnDir,
                getGeneNames=FALSE, 
                pop1alias="mel", pop2alias="sim")


##### all Courtney new alignments, before removing seqs with Ns, unpolarized tests
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

## unpolarized test on all Arps alignments - works
Arps_inclNseqs_results <- lapply (names(Arps_inclNseqs_alnFiles), function(x){
    doMKtest(Arps_inclNseqs_alnFiles[x], 
          pop1seqs=Arps_inclNseqs_populationNames[[x]][["pop1"]], pop1alias="mel",
          pop2seqs=Arps_inclNseqs_populationNames[[x]][["pop2"]], pop2alias="sim" )
})
names(Arps_inclNseqs_results) <- names(Arps_inclNseqs_alnFiles)

## combine results from all Arps alignments
Arps_inclNseqs_results_all <- combineMKresults(Arps_inclNseqs_results, 
    outFile="Arps_inclNseqs_allMKresults_conservative.xlsx",
    outDir=Arps_inclNseqs_alnDir,
    getGeneNames=FALSE, 
    pop1alias="mel", pop2alias="sim")




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

## unpolarized test on all Arps alignments - works
Arps_noNseqs_results <- lapply (names(Arps_noNseqs_alnFiles), function(x){
    doMKtest(Arps_noNseqs_alnFiles[x], 
          pop1seqs=Arps_noNseqs_populationNames[[x]][["pop1"]], pop1alias="mel",
          pop2seqs=Arps_noNseqs_populationNames[[x]][["pop2"]], pop2alias="sim" )
})
names(Arps_noNseqs_results) <- names(Arps_noNseqs_alnFiles)

## combine results from all Arps alignments
Arps_noNseqs_results_all <- combineMKresults(Arps_noNseqs_results, 
    outFile="Arps_noNseqs_allMKresults_conservative.xlsx",
    outDir=Arps_noNseqs_alnDir,
    getGeneNames=FALSE, 
    pop1alias="mel", pop2alias="sim")


##### add the corresponding Arp53D results to those for the other Arps

ArpsPlusArp53D_noNseqs_results <- c( Arps_noNseqs_results, 
                                     Arp53D_allResults["noN_justZI"])
names(ArpsPlusArp53D_noNseqs_results) <- gsub("noN_justZI","Arp53D",names(ArpsPlusArp53D_noNseqs_results))

ArpsPlusArp53D_noNseqs_results_all <- combineMKresults(ArpsPlusArp53D_noNseqs_results, 
                                    outFile="Arps_plus_Arp53D_allMKresults_noNseqs.xlsx",
                                    outDir="data/Arp53D",
                                    getGeneNames=FALSE, 
                                    pop1alias="mel", pop2alias="sim")


ArpsPlusArp53D_inclNseqs_results <- c( Arps_inclNseqs_results, Arp53D_allResults["inclN_justZI"])
names(ArpsPlusArp53D_inclNseqs_results) <- gsub("inclN_justZI","Arp53D",names(ArpsPlusArp53D_inclNseqs_results))
ArpsPlusArp53D_inclNseqs_results_all <- combineMKresults(ArpsPlusArp53D_inclNseqs_results, 
                                    outFile="Arps_plus_Arp53D_allMKresults_inclNseqs.xlsx",
                                    outDir="data/Arp53D",
                                    getGeneNames=FALSE, 
                                    pop1alias="mel", pop2alias="sim")





##### Gustavo's alignments, unpolarized tests

gustavo_alnDir <- "data/Gustavo/originalFilesAsFasta"

gustavo_alnFiles <- list.files(gustavo_alnDir, pattern=".fa$", full.names=TRUE)
names(gustavo_alnFiles) <- gsub(paste(gustavo_alnDir, "/", sep=""), "", gustavo_alnFiles)
names(gustavo_alnFiles) <- gsub(".fa$", "", names(gustavo_alnFiles))

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
    populationNames[["out"]][[1]] <- grep("outgroup", names(x), value=TRUE)
    return(populationNames)
})
# sapply(gustavo_populationNames, function(x) { sapply(x, length)})

## unpolarized test on all Gustavo alignments - works
gustavo_results <- lapply (names(gustavo_alnFiles), function(x){
    doMKtest(gustavo_alnFiles[x], 
             outDir="data/Gustavo/MKresults_unpolarized",
             pop1seqs=gustavo_populationNames[[x]][["pop1"]], pop1alias="buz",
             pop2seqs=gustavo_populationNames[[x]][["pop2"]], pop2alias="sma" )
})
names(gustavo_results) <- names(gustavo_alnFiles)

## combine results from all Gustavo alignments
gustavo_results_all <- combineMKresults(gustavo_results, 
                                        outFile="gustavo_allMKresults_unpolarized.xlsx",
                                        outDir="data/Gustavo/MKresults_unpolarized",
                                        getGeneNames=FALSE, 
                                        pop1alias="buz", pop2alias="sma")

## polarized test on all Gustavo alignments - works
gustavo_results_polarized <- lapply (names(gustavo_alnFiles), function(x){
    doMKtest(gustavo_alnFiles[x], 
             outDir="data/Gustavo/MKresults_polarized",
             pop1seqs=gustavo_populationNames[[x]][["pop1"]], pop1alias="buz",
             pop2seqs=gustavo_populationNames[[x]][["pop2"]], pop2alias="sma",
             outgroupSeqs=gustavo_populationNames[[x]][["out"]],
             polarize=TRUE )
})
names(gustavo_results_polarized) <- names(gustavo_alnFiles)

## combine results from all Gustavo alignments
gustavo_results_polarized_all <- combineMKresults(gustavo_results_polarized, 
                                        outFile="gustavo_allMKresults_polarized.xlsx",
                                        outDir="data/Gustavo/MKresults_polarized",
                                        getGeneNames=FALSE, 
                                        pop1alias="buz", pop2alias="sma")


pdf(height=5,width=15,file="data/Gustavo/MKresults_polarized/gustavo_allMKresults_polarized.pdf")
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

# a single plot
#temp <- plotMKpositions(gustavo_results_polarized[["Cid6pop-onlyHFD"]][["positions"]], 
#                        plotPolarizedPop1=TRUE, plotPolarizedPop2=TRUE, 
#                        title="Cid6pop-onlyHFD", pop1alias="buz", pop2alias="sma")



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



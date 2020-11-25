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
MKresults_websiteExample <- doMKtest(
          "test_data/MKTwebsiteExample/MKTwebsite_testAln/MKTwebsite_testAln.fa", 
          pop1seqs=c("pongo1","pongo2"), pop1alias="pongo",
          pop2seqs=c("trachy1","trachy2"), pop2alias="trachy",
          writeAncFasta=TRUE, writeMKoutput=TRUE)

### and make a plot of where the changes are
pdf(height=5,width=15,file="test_data/MKTwebsiteExample/MKTwebsite_testAln/MKTwebsite_testAln_MKplot.pdf")
plotMKpositions(MKresults_websiteExample[["positions"]], 
                title="MKresults websiteExample", pop1alias="pongo", pop2alias="trachy", 
                setNumPlots=FALSE)
dev.off()

######## test using same alignment as MK website, but now with seqs including large stretches of Ns ior gaps added.  On the website, this means NO changes are counted, but I still count those positions

MKresults_websiteExample_addNseqs <- doMKtest(
    "test_data/MKTwebsiteExample/MKTwebsite_testAln_addNseqs/MKTwebsite_testAln_addNseqs.fa", 
    pop1seqs=c("pongo1","pongo2","pongo3","pongo4"), pop1alias="pongo",
    pop2seqs=c("trachy1","trachy2"), pop2alias="trachy",
    writeAncFasta=FALSE, writeMKoutput=TRUE)

MKresults_websiteExample_addGapseqs <- doMKtest(
    "test_data/MKTwebsiteExample/MKTwebsite_testAln_addNseqs/MKTwebsite_testAln_addGapSeqs.fa", 
    pop1seqs=c("pongo1","pongo2","pongo3","pongo4"), pop1alias="pongo",
    pop2seqs=c("trachy1","trachy2"), pop2alias="trachy",
    writeAncFasta=FALSE, writeMKoutput=TRUE)

temp_aln <- readBStringSet("test_data/MKTwebsiteExample/MKTwebsite_testAln_addNseqs/MKTwebsite_testAln_addNseqs.fa")
names(temp_aln) <- sapply(strsplit(names(temp_aln)," "),"[[",1)
temp_aln_df <- as.data.frame(as.matrix(temp_aln), stringsAsFactors=FALSE)
temp_aln_df_freqs <- getACGTfreqs(tabulateDF(temp_aln_df))



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



## combine results from all test3 alignments
test3_results_all <- combineMKresults(test3_results, 
                                       outFile="test3_results_freqFilter_20pc.xlsx",
                                       outDir=test3_dir,
                                       getGeneNames=FALSE)



#test3_results[["summary"]]
#head(test3_results[["positions"]])




### https://en.wikipedia.org/wiki/McDonald%E2%80%93Kreitman_test
# alpha represents the proportion of substitutions driven by positive selection.
# neutrality index: >1 indicates negative selection; <1 indicates positive selection

### set working directory
# mac location
#setwd("~/Desktop/mac_workStuff/mac_MKtests/MKtests_JY/")
# mac location - server
setwd("/Volumes/malik_h/user/jayoung/MKtest/MKtests_JY")
# fast location
#setwd("/fh/fast/malik_h/user/jayoung/MKtest/forGithub/")

### load functions
source("scripts/MKfunctions.R")

######## test using same alignment as MK website:
MKresults_websiteExample_outDir <- "test_data/MKTwebsiteExample/MKTwebsite_testAln/MKTwebsite_testAln_testOutput"

MKresults_websiteExample <- doMKtest(
          myAlnFile="test_data/MKTwebsiteExample/MKTwebsite_testAln/MKTwebsite_testAln.fa", 
          outDir=MKresults_websiteExample_outDir, 
          pop1seqs=c("pongo1","pongo2"), pop1alias="pongo",
          pop2seqs=c("trachy1","trachy2"), pop2alias="trachy",
          writeAncFasta=TRUE, writeMKoutput=TRUE)

### and make a plot of where the changes are
pdf(height=5,width=15,
    file=paste(MKresults_websiteExample_outDir,"MKTwebsite_testAln_MKplot.pdf",sep="/"))
plotMKpositions(MKresults_websiteExample[["positions"]], 
                title="MKresults websiteExample", pop1alias="pongo", pop2alias="trachy", 
                setNumPlots=FALSE)
dev.off()

######## test using same alignment as MK website, but now with seqs including large stretches of Ns added.  On the website, this means NO changes are counted, but I still count those positions

MKresults_websiteExample_addNseqs_outDir <- "test_data/MKTwebsiteExample/MKTwebsite_testAln_addNseqs/MKTwebsite_testAln_addNseqs_testOutput"

MKresults_websiteExample_addNseqs <- doMKtest(
    "test_data/MKTwebsiteExample/MKTwebsite_testAln_addNseqs/MKTwebsite_testAln_addNseqs.fa", 
    outDir=MKresults_websiteExample_addNseqs_outDir, 
    pop1seqs=c("pongo1","pongo2","pongo3","pongo4"), pop1alias="pongo",
    pop2seqs=c("trachy1","trachy2"), pop2alias="trachy",
    writeAncFasta=FALSE, writeMKoutput=TRUE)


# read in that alignment and take a look at stats
temp_aln <- readBStringSet("test_data/MKTwebsiteExample/MKTwebsite_testAln_addNseqs/MKTwebsite_testAln_addNseqs.fa")
names(temp_aln) <- sapply(strsplit(names(temp_aln)," "),"[[",1)
temp_aln_df <- as.data.frame(as.matrix(temp_aln), stringsAsFactors=FALSE)
temp_aln_df_freqs <- getACGTfreqs(tabulateDF(temp_aln_df))


######## test using same alignment as MK website, but now with seqs including large stretches of gaps added.  On the website, this means NO changes are counted, but I still count those positions

MKresults_websiteExample_addGapseqs_outDir <- "test_data/MKTwebsiteExample/MKTwebsite_testAln_addNseqs/MKTwebsite_testAln_addGapSeqs_testOutput"

MKresults_websiteExample_addGapseqs <- doMKtest(
    "test_data/MKTwebsiteExample/MKTwebsite_testAln_addNseqs/MKTwebsite_testAln_addGapSeqs.fa", 
    outDir=MKresults_websiteExample_addGapseqs_outDir, 
    pop1seqs=c("pongo1","pongo2","pongo3","pongo4"), pop1alias="pongo",
    pop2seqs=c("trachy1","trachy2"), pop2alias="trachy",
    writeAncFasta=FALSE, writeMKoutput=TRUE)



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


####### test4 - a weird polymorphism issue that came up in one of Ching-Ho's alignments (BigH1_CDS_with_sim.fasta)
# pop1_seq1   GCG
# pop1_seq2   TCG
# pop1_seq3   GAG
# pop2_seq1   GAG


# In this case, the ancestral codon of pop1 is most likely to be GAG.
# The A->C polymorphism at position 2 is counted OK as a non-synonymous polymorphism
# The G->T polymorphism at position 1 created a problem for one version of my code - if it occurred on the background of the ancestral codon it would create a stop codon (GAG->TAG), which I don't know how to count So it must have occurred on the 'intermediate' codon (GCG->TCG), and therefore is also a non-synonymous codon.
# I have fixed the initial problem, where it simply broke the code. However, I don't think I am counting it right now. The MKT website counts this as 2 non-synonymous polymorphisms. 
# I count changes from before codon -> after codon
# I'm using the inferred ancestral codon as before. But perhaps I need to count over all possible before codons for polymorphisms.  All pairwise combinations of codons ??  see if the MKT website says how they do it



test4_alnFile <- "test_data/test_unusualPolymorphism/test_unusualPolymorphism.fa"


test4_result <- doMKtest(myAlnFile=test4_alnFile, 
                         outfileStem="test4", 
                         pop1seqs=c("pop1_seq1","pop1_seq2","pop1_seq3"),
                         pop2seqs=c("pop2_seq1"), 
                         quiet=TRUE )


test4_result[["positions"]][,c("pop1_anc","pop2_anc",
                               "pop1_vs_pop2_Dn", "pop1_vs_pop2_Ds", 
                               "pop1_Pn", "pop1_Ps" )]

test4_result[["positions"]][,c("pop1_A", "pop1_C", "pop1_G", "pop1_T", 
                               "pop2_A", "pop2_C", "pop2_G", "pop2_T" )]


test4_aln <- readBStringSet(test4_alnFile)


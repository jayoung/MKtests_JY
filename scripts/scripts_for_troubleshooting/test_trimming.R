source("scripts/MKfunctions.R")

######## test using same alignment as MK website:

alnFile <- "test_data/test_trimming/MKTwebsite_testAln.fa"

aln <- readBStringSet(alnFile)

MKresults_full <- doMKtest(
    myAlnFile=alnFile, 
    pop1seqs=c("pongo1","pongo2"), pop1alias="pongo",
    pop2seqs=c("trachy1","trachy2"), pop2alias="trachy")

MKresults_codon21 <- doMKtest(
    myAlnFile=alnFile, 
    pop1seqs=c("pongo1","pongo2"), pop1alias="pongo",
    pop2seqs=c("trachy1","trachy2"), pop2alias="trachy",
    regionStartAA=21, regionEndAA=21)

MKresults_chunk1 <- doMKtest(
    myAlnFile=alnFile, 
    pop1seqs=c("pongo1","pongo2"), pop1alias="pongo",
    pop2seqs=c("trachy1","trachy2"), pop2alias="trachy",
    regionEndAA=21)


MKresults_chunk2 <- doMKtest(
    myAlnFile=alnFile, 
    pop1seqs=c("pongo1","pongo2"), pop1alias="pongo",
    pop2seqs=c("trachy1","trachy2"), pop2alias="trachy",
    regionStartAA=22)

MKresults_combined <- combineMKresults(list(MKresults_full=MKresults_full,
                                            MKresults_chunk1=MKresults_chunk1, 
                                            MKresults_chunk2=MKresults_chunk2,
                                            MKresults_codon21=MKresults_codon21),
                                       outFile="test_trimming.xlsx",
                                       outDir="test_data/test_trimming")

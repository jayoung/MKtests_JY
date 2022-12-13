### load Janet's MK functions
source("scripts/MKfunctions.R")


### define an input alignment file:
alnFile <- "test_data/MKTwebsiteExample/MKTwebsite_testAln.fa"

### read in the alignment file, so that we can figure out the names of the sequences in each population
aln <- readBStringSet(alnFile)
# fix the seqnames in the alignment - only want the bits before the space:
names(aln) <- sapply( strsplit(names(aln), " "), "[[", 1)

### figure out the names in each population. This step will be different for your alignment.
pongo_names <- grep("pongo", names(aln), value=TRUE)
trachy_names <- grep("trachy", names(aln), value=TRUE)

### do the MK test - simplest example
MKresults_websiteExample <- doMKtest(
          myAlnFile=alnFile, 
          pop1seqs=pongo_names, 
          pop2seqs=trachy_names)

# the output of the doMKtest function is a list, comprising two tables.
#    list element 1 ("summary")  contains the MK test results
#    list element 2 ("positions")   contains a per-position analysis of the alignment. This helps understand what's going on with the alignment, and/or with troubleshooting



### do the same MK test and save an Excel file of the results (we add the writeMKoutput, outDir and outfileStem options)
# (if we also use aliases, as below, it'll help us interpret the column names in the output)
MKresults_websiteExample <- doMKtest(
    myAlnFile=alnFile, 
    writeMKoutput=TRUE, outDir = "test_data/MKTwebsiteExample", outfileStem = "MKresults",
    pop1seqs=pongo_names, pop1alias = "pongo",
    pop2seqs=trachy_names, pop2alias = "trachy")
# as well as the list in R, we now also have an Excel output file called "MKTwebsite_testAln.MK.xlsx". There are two tabs, representing the two data.frames in the list output


### if we've performed more than one MK test, we can combine the results like this:
MKresults_combined <- combineMKresults(list(results1=MKresults_websiteExample,
                                            results2=anotherMKresult),
                                       outFile="combinedMKresults.xlsx",
                                       outDir=".")
## the R object that's returned will be a single data.frame, combining the results.
# the Excel file will have that table summarizing all the MK results, and one additional tab for each MK alignment showing the positions table

### we can make a plot showing where the fixed and polymorphic changes are:
pdf(height=5,width=15, file="test_data/MKTwebsiteExample/MKplot.pdf")
plotMKpositions(MKresults_websiteExample[["positions"]], 
                title="MKresults websiteExample", pop1alias="pongo", pop2alias="trachy", 
                setNumPlots=FALSE)
dev.off()



#### useful options 

## filter rare variants: use the filterRareAlleles/alleleFreqThreshold options
# (the example alignment only has two individuals for each population, so the filtering makes no sense for this example)
MKresults_websiteExample_removeRareVariants <- doMKtest(
    myAlnFile=alnFile, 
    pop1seqs=pongo_names, pop1alias = "pongo",
    pop2seqs=trachy_names, pop2alias = "trachy",
    filterRareAlleles=TRUE, alleleFreqThreshold=0.05)

## only look at part of an alignment: use the regionStartAA/regionEndAA options to specify the start/end positions (in amino acids) of the region we want to look at
MKresults_websiteExample_region1 <- doMKtest(
    myAlnFile=alnFile, 
    pop1seqs=pongo_names, pop1alias = "pongo",
    pop2seqs=trachy_names, pop2alias = "trachy",
    regionStartAA=1, regionEndAA=21)


##### polarized MK tests

## add polarize=TRUE and supply outgroup seqnames via the outgroupSeqs option

## example data - CG17802 (Nicknack) in D.melanogaster and D.simulans with D. yakuba outgroup
# see https://elifesciences.org/articles/63368

nnk_alnFile <- "test_data/Nnk_Kasinathan_2020/802.plusOutgroup_aln1_NT.edit.noStops.degapcodon.fas"

nnk_aln <- readBStringSet(nnk_alnFile)

nnk_Dmel_names <- names(nnk_aln)[3:length(nnk_aln)]
nnk_Dsim_name <- names(nnk_aln)[2]
nnk_Dyak_names <- names(nnk_aln)[1]

nnk_MKresults_polarized <- doMKtest(
    myAlnFile=nnk_alnFile, 
    pop1seqs=nnk_Dmel_names, pop1alias = "Dmel",
    pop2seqs=nnk_Dsim_name, pop2alias = "Dsim",
    polarize=TRUE,
    outgroupSeqs=nnk_Dyak_names )






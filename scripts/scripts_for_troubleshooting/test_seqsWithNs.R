## I think seqs with Ns mess up the plot, for polymorphisms (but not fixed changes).
## troubleshoot that.

library(here)

source("scripts/MKfunctions.R")

## define input files
alnFiles <- list.files(here("test_data/test_seqsWithNs"), 
                       pattern=".fa$", full.names=TRUE)
names(alnFiles) <- c("withNs", "withoutNs","singleCodon_withNs", "singleCodon_withoutNs")

## read alignments and sort out names
alns <- lapply(alnFiles, function(x) {
    y <- readBStringSet(x)
    names(y) <- sapply(strsplit(names(y), " "), "[[", 1)
    return(y)
})

## get IDs of individuals in each population
pops <- lapply(alns, function(x) {
    species <- gsub("\\d","",names(x))
    populations <- split(names(x), species)
    return(populations)
})


#### single codon aln, without Ns
MKresults_withoutNs <- doMKtest(alnFiles[["singleCodon_withoutNs"]], 
                                pop1seqs=pops[["singleCodon_withoutNs"]][["pongo"]],
                                pop2seqs=pops[["singleCodon_withoutNs"]][["trachy"]])

#### single codon aln, with Ns
MKresults_withNs <- doMKtest(alnFiles[["singleCodon_withNs"]], 
                             pop1seqs=pops[["singleCodon_withNs"]][["pongo"]],
                             pop2seqs=pops[["singleCodon_withNs"]][["trachy"]])


#### example alignment from MKT website
MKresults_withoutNs <- doMKtest(alnFiles[["withoutNs"]], 
                                pop1seqs=pops[["withoutNs"]][["pongo"]],
                                pop2seqs=pops[["withoutNs"]][["trachy"]])

#### example alignment from MKT website, with Ns
MKresults_withNs <- doMKtest(alnFiles[["withNs"]], 
                             pop1seqs=pops[["withNs"]][["pongo"]],
                             pop2seqs=pops[["withNs"]][["trachy"]])



### and make a plot of where the changes are
pdf(height=5,width=15,
    file=here("test_data/test_seqsWithNs/MKresults_withoutNs_MKplot.pdf"))

plotMKpositions(MKresults_withoutNs[["positions"]], 
                title="MKresults websiteExample without Ns", pop1alias="pongo", pop2alias="trachy", 
                setNumPlots=FALSE)
dev.off()




### and make a plot of where the changes are
pdf(height=5,width=15,
    file=here("test_data/test_seqsWithNs/MKresults_withNs_MKplot.pdf"))

plotMKpositions(MKresults_withNs[["positions"]], 
                title="MKresults websiteExample with Ns", pop1alias="pongo", pop2alias="trachy", 
                setNumPlots=FALSE)
dev.off()


### check with Risa's original alignment
## define input files
risaNalnFile <- here("data/Risa/2022_Nov29/Janet_abo_alignment.fasta")
risaNaln <- readBStringSet(risaNalnFile)

risaNalnPops <- list()
risaNalnPops[["mel"]] <- names(risaNaln)[1:1102]
risaNalnPops[["sim"]] <- names(risaNaln)[1103]


#### example alignment from MKT website, with Ns
risaNalnPops_MKresults <- doMKtest(risaNalnFile, 
                                   pop1seqs=risaNalnPops[["mel"]],
                                   pop2seqs=risaNalnPops[["sim"]])


### and make a plot of where the changes are
pdf(height=5,width=15,
    file=here("test_data/test_seqsWithNs/risaNalnPops_MKresults_MKplot.pdf"))

plotMKpositions(risaNalnPops_MKresults[["positions"]], 
                title="Risa mel-simabo alignment with Ns", pop1alias="mel", pop2alias="sim", 
                setNumPlots=FALSE)
dev.off()



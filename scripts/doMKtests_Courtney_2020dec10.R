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


## before using R I did this:
# cd /Users/jayoung/Desktop/mac_workStuff/mac_MKtests/MKtests_JY/data/Courtney_2020_Dec10/melAndSim
# cat ../fromCourtney/dhc16F/Dsim_dhc16F_CDS.fasta  ../fromCourtney/dhc16F/Nucleotide\ alignment_minus\ introns_minus\ gene\ model.fasta  > dhc16F_aln.fa
# cat ../fromCourtney/klp3a/Dsim_klp3a_XM_016182791.1.fasta ../fromCourtney/klp3a/Nucleotide\ alignment\ 2_minus\ introns_klp3a_minus_gene_model.fasta > klp3a_aln.fa
# cat ../fromCourtney/lamin_b/dsim_XM_016184192.1_laminb.fasta  ../fromCourtney/lamin_b/Nucleotide\ alignment\ 2_minus_introns_minus_gene_model.fasta > lamin_b_aln.fa
#     edit each alignment by hand, and remove the dmel reference seq that's not part of the ZI population


alnDir <- "data/Courtney_2020_Dec10/melAndSim"

alnFiles <- list.files(alnDir, full.names=TRUE, pattern=".fa")
names(alnFiles) <- gsub("_aln.fa","",sapply(strsplit(alnFiles, "/"), "[[", 4))


alns <- lapply(alnFiles, function(x) {
    aln <- readBStringSet(x)
    names(aln) <- sapply(strsplit(names(aln)," "),"[[",1)
    return(aln)
})

# get names of seqs within each population
populationNames <- lapply(alns, function(x){
    populationNames <- list()
    populationNames[["pop1"]] <- grep("ZI", names(x), value=TRUE)
    populationNames[["pop2"]] <- grep("sim", names(x), value=TRUE)
    return(populationNames)
})
sapply(populationNames, function(x) { sapply(x, length)})

# get names of seqs within each population that do not contain any Ns
populationNames_noNs <- lapply(alns, function(x){
    ## remove seqs containing Ns
    alnFilt <- x[ which(!grepl("N", as.character(x), ignore.case=TRUE)) ]
    populationNames_noNs <- list()
    populationNames_noNs[["pop1"]] <- grep("ZI", names(alnFilt), value=TRUE)
    populationNames_noNs[["pop2"]] <- grep("sim", names(alnFilt), value=TRUE)
    return(populationNames_noNs)
})
sapply(populationNames_noNs, function(x) { sapply(x, length)})



### do various flavors of MK test on those three alignments.  all will be unpolarized, as there's no outgroup
# all seqs, all SNPs
# all seqs, remove alleles <= 5% frequency
# remove seqs containing Ns, all alleles
# remove seqs containing Ns, remove alleles <= 5% frequency

# unpolarized all seqs, all SNPs
results_allSeqsAllAlleles <- lapply (names(alns), function(x){
    doMKtest(alnFiles[x], 
             pop1seqs=populationNames[[x]][["pop1"]], pop1alias="mel",
             pop2seqs=populationNames[[x]][["pop2"]], pop2alias="sim" )
})
names(results_allSeqsAllAlleles) <- names(alns)

results_allSeqsAllAlleles_combined <- combineMKresults(results_allSeqsAllAlleles, 
                    outFile="threeMoreGenes_allSeqs_allAlleles.xlsx",
                    outDir=alnDir,
                    getGeneNames=FALSE, 
                    pop1alias="mel", pop2alias="sim")

# unpolarized all seqs, remove alleles <= 5% frequency
results_allSeqs_freqFilt <- lapply (names(alns), function(x){
    doMKtest(alnFiles[x], 
             pop1seqs=populationNames[[x]][["pop1"]], pop1alias="mel",
             pop2seqs=populationNames[[x]][["pop2"]], pop2alias="sim",
             filterRareAlleles=TRUE, alleleFreqThreshold=0.05 )
})
names(results_allSeqs_freqFilt) <- names(alns)

results_allSeqs_freqFilt_combined <- combineMKresults(results_allSeqs_freqFilt, 
                    outFile="threeMoreGenes_allSeqs_freqFilt.xlsx",
                    outDir=alnDir,
                    getGeneNames=FALSE, 
                    pop1alias="mel", pop2alias="sim")

# remove seqs containing Ns, all alleles
results_filtSeqs_allAlleles <- lapply (names(alns), function(x){
    doMKtest(alnFiles[x], 
             pop1seqs=populationNames_noNs[[x]][["pop1"]], pop1alias="mel",
             pop2seqs=populationNames_noNs[[x]][["pop2"]], pop2alias="sim" )
})
names(results_filtSeqs_allAlleles) <- names(alns)

results_filtSeqs_allAlleles_combined <- combineMKresults(results_filtSeqs_allAlleles, 
                    outFile="threeMoreGenes_filtSeqs_allAlleles.xlsx",
                    outDir=alnDir,
                    getGeneNames=FALSE, 
                    pop1alias="mel", pop2alias="sim")

# remove seqs containing Ns, remove alleles <= 5% frequency
results_filtSeqs_filtAlleles <- lapply (names(alns), function(x){
    doMKtest(alnFiles[x], 
             pop1seqs=populationNames_noNs[[x]][["pop1"]], pop1alias="mel",
             pop2seqs=populationNames_noNs[[x]][["pop2"]], pop2alias="sim",
             filterRareAlleles=TRUE, alleleFreqThreshold=0.05  )
})
names(results_filtSeqs_filtAlleles) <- names(alns)

results_filtSeqs_filtAlleles_combined <- combineMKresults(results_filtSeqs_filtAlleles, 
                    outFile="threeMoreGenes_filtSeqs_freqFilt.xlsx",
                    outDir=alnDir,
                    getGeneNames=FALSE, 
                    pop1alias="mel", pop2alias="sim")


#### plots

plotWidth <- 22
plotHeight <- 6

### all seqs, without allele freq filter
pdf(height=plotHeight,width=plotWidth, 
    file="data/Courtney_2020_Dec10/melAndSim/threeMoreGenes_allSeqs_allAlleles_MKplots.pdf")
par(mfrow=c(3,1), oma=c(0,0,2,0))
for (x in names(results_allSeqsAllAlleles)) {
    cat("plotting",x,"\n")
    plotMKpositions(results_allSeqsAllAlleles[[x]][["positions"]], 
                    title=x, 
                    setNumPlots=FALSE)
}
title(main="all seqs (including seqs with Ns), no allele frequency filter", 
      outer=TRUE, cex.main=2)
dev.off()

### all seqs, 5% allele freq filter
pdf(height=plotHeight,width=plotWidth, 
    file="data/Courtney_2020_Dec10/melAndSim/threeMoreGenes_allSeqs_freqFilt_MKplots.pdf")
par(mfrow=c(3,1), oma=c(0,0,2,0))
for (x in names(results_allSeqs_freqFilt)) {
    cat("plotting",x,"\n")
    plotMKpositions(results_allSeqs_freqFilt[[x]][["positions"]], 
                    title=x, 
                    setNumPlots=FALSE)
}
title(main="all seqs (including seqs with Ns), remove alleles with freq <= 5%", 
      outer=TRUE, cex.main=2)
dev.off()

### only seqs without Ns, without allele freq filter
pdf(height=plotHeight,width=plotWidth, 
    file="data/Courtney_2020_Dec10/melAndSim/threeMoreGenes_filtSeqs_allAlleles_MKplots.pdf")
par(mfrow=c(3,1), oma=c(0,0,2,0))
for (x in names(results_filtSeqs_allAlleles)) {
    cat("plotting",x,"\n")
    plotMKpositions(results_filtSeqs_allAlleles[[x]][["positions"]], 
                    title=x, 
                    setNumPlots=FALSE)
}
title(main="only seqs without Ns, no allele frequency filter", 
      outer=TRUE, cex.main=2)
dev.off()

### only seqs without Ns, 5% allele freq filter
pdf(height=plotHeight,width=plotWidth, 
    file="data/Courtney_2020_Dec10/melAndSim/threeMoreGenes_filtSeqs_freqFilt_MKplots.pdf")
par(mfrow=c(3,1), oma=c(0,0,2,0))
for (x in names(results_filtSeqs_filtAlleles)) {
    cat("plotting",x,"\n")
    plotMKpositions(results_filtSeqs_filtAlleles[[x]][["positions"]], 
                    title=x, 
                    setNumPlots=FALSE)
}
title(main="only seqs without Ns, remove alleles with freq <= 5%", 
      outer=TRUE, cex.main=2)
dev.off()



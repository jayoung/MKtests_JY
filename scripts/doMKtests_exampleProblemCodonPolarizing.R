## figuring out a bug.   Risa's mel-sim polarized test revealed a potential problem.
library(Biostrings)
library(here)

### load functions
source(here("scripts/MKfunctions.R"))


#### example codon change, that was being treated oddly by my script before Dec 7 2022

# Risa and I noticed this when digging in to her Abo results

# it's a synonymous change, but when we try to polarize, we cannot tell which branch it occurred on. Should not be a count on either branch. In the pre-Dec 7 code it always gets assigned to pop2, and therefore depended on the order in which I specify the two populations. Turns out that's to do with how I select the counts when there's >1 possible before-after combination of codons. I WAS choosing NS,S counts that minimized NS, and if there was a tie I arbitrarily chose the first set of counts.  Now I choose to ALSO minimize S, if there's a tie on NS counts.


## set up a very small alignment - just one codon
aln <- character() 
aln["popA_ind1"] <- "ACA"
aln["popA_ind2"] <- "ACA"
aln["popB_ind1"] <- "ACC"
aln["popB_ind2"] <- "ACC"
aln["outgroup"]  <- "ACT"
aln <- DNAStringSet(aln)

## save as fasta
alnFile <- here("test_data/test_ambiguousCodon/test_ambiguousCodon.fa") 
writeXStringSet(aln, alnFile)

## set up the population names
populationIDs <- list()
populationIDs[["popA"]] <- names(aln)[1:2]
populationIDs[["popB"]] <- names(aln)[3:4]
populationIDs[["out"]] <- names(aln)[5]

MKresults <- list()

MKresults[["unpolarizedA"]] <- doMKtest(alnFile,           
                                        #outDir=here("test_data/test_ambiguousCodon/unpolarizedA"), 
                                        pop1seqs=populationIDs[["popA"]], pop1alias="popA",
                                        pop2seqs=populationIDs[["popB"]], pop2alias="popB"#,
                                        #writeMKoutput=TRUE
)

MKresults[["unpolarizedB"]] <- doMKtest(alnFile, 
                                        pop1seqs=populationIDs[["popB"]], pop2alias="popB",
                                        pop2seqs=populationIDs[["popA"]], pop1alias="popA")


MKresults[["polarizedA"]] <- doMKtest(alnFile, 
                                      pop1seqs=populationIDs[["popA"]], pop1alias="popA",
                                      pop2seqs=populationIDs[["popB"]], pop2alias="popB",
                                      polarize=TRUE, 
                                      outgroupSeqs = list(populationIDs[["out"]]),
                                      extraVerbose=FALSE)


MKresults[["polarizedB"]] <- doMKtest(alnFile, 
                                      pop1seqs=populationIDs[["popB"]], pop2alias="popB",
                                      pop2seqs=populationIDs[["popA"]], pop1alias="popA",
                                      polarize=TRUE, 
                                      outgroupSeqs = list(populationIDs[["out"]]),
                                      extraVerbose=FALSE )

## combine and save to Excel files
MKresults_unpolarized <- combineMKresults(MKresults[1:2], 
                                          outFile="test_ambiguousCodon_MKresults_unpolarized.xlsx",
                                          outDir=here("test_data/test_ambiguousCodon"),
                                          getGeneNames=FALSE)

MKresults_polarized <- combineMKresults(MKresults[3:4], 
                                        outFile="test_ambiguousCodon_MKresults_polarized.xlsx",
                                        outDir=here("test_data/test_ambiguousCodon"),
                                        getGeneNames=FALSE,
                                        extraVerbose=TRUE)



#### Risa's alignments:

# her older results are in the oldResults_risa_beforeFixConservative folder

## Risa said:
# For the “melpop_simpop” file, seqs 1-36 are sim, seqs 37-650 are mel, and seq 651 is yak.
# 
# The “melpop_simpop” files contain 2 more sim sequences than the “simpop” file. These are 1) the reference sim sequence I found on NCBI, and 2) the sim sequence that Ching-Ho used for the original MK test (which is not included in the new sim population data he gave me).

## in popDat_bothSpecies, but NOT popDat_justSim, these positions were classified as fixed changes in the sim branch:   151, 208, 373, 658, 1603

## in popDat_justSim, but NOT popDat_bothSpecies, these positions were classified as fixed changes in the sim branch:   949, 1588

risaFiles <- list()
risaFiles[["popDat_justSim"]] <- here("data/Risa/2022_Dec7/abo_simpop_mel_yak.fasta")
risaFiles[["popDat_bothSpecies"]] <- here("data/Risa/2022_Dec7/abo_melpop_simpop_yak.fasta")

risa_alns <- lapply(risaFiles, readDNAStringSet)

## set up population names
risa_popIDs <- list()

risa_popIDs[["popDat_justSim"]] <- list()
risa_popIDs[["popDat_justSim"]][["sim"]] <- grep("CM015606", names(risa_alns[["popDat_justSim"]]), value=TRUE)
risa_popIDs[["popDat_justSim"]][["mel"]] <- grep("melanogaster", names(risa_alns[["popDat_justSim"]]), value=TRUE)
risa_popIDs[["popDat_justSim"]][["yak"]] <- grep("yakuba", names(risa_alns[["popDat_justSim"]]), value=TRUE)

risa_popIDs[["popDat_bothSpecies"]] <- list()
risa_popIDs[["popDat_bothSpecies"]][["sim"]] <- grep("CM015606|sim", names(risa_alns[["popDat_bothSpecies"]]), value=TRUE)
risa_popIDs[["popDat_bothSpecies"]][["yak"]] <- grep("yakuba", names(risa_alns[["popDat_bothSpecies"]]), value=TRUE)
## mel is everything else:
risa_popIDs[["popDat_bothSpecies"]][["mel"]] <- setdiff(
    names(risa_alns[["popDat_bothSpecies"]]), 
    c(risa_popIDs[["popDat_bothSpecies"]][["sim"]],
      risa_popIDs[["popDat_bothSpecies"]][["yak"]])) 
## sanity check
# lapply(risa_popIDs, function(x) { sapply(x, length)})

## MK tests
risa_MKresults <- list()

## popDat_justSim
risa_MKresults[["popDat_justSim_unpolarizedA"]] <- doMKtest(
    risaFiles[["popDat_justSim"]], 
    pop1seqs=risa_popIDs[["popDat_justSim"]][["sim"]], pop1alias="sim",
    pop2seqs=risa_popIDs[["popDat_justSim"]][["mel"]], pop2alias="mel")

risa_MKresults[["popDat_justSim_unpolarizedB"]] <- doMKtest(
    risaFiles[["popDat_justSim"]], 
    pop1seqs=risa_popIDs[["popDat_justSim"]][["mel"]], pop1alias="mel",
    pop2seqs=risa_popIDs[["popDat_justSim"]][["sim"]], pop2alias="sim")

risa_MKresults[["popDat_justSim_polarizedA"]] <- doMKtest(
    risaFiles[["popDat_justSim"]], 
    pop1seqs=risa_popIDs[["popDat_justSim"]][["sim"]], pop1alias="sim",
    pop2seqs=risa_popIDs[["popDat_justSim"]][["mel"]], pop2alias="mel",
    polarize=TRUE, 
    outgroupSeqs = list(risa_popIDs[["popDat_justSim"]][["yak"]]))

risa_MKresults[["popDat_justSim_polarizedB"]] <- doMKtest(
    risaFiles[["popDat_justSim"]], 
    pop1seqs=risa_popIDs[["popDat_justSim"]][["mel"]], pop1alias="mel",
    pop2seqs=risa_popIDs[["popDat_justSim"]][["sim"]], pop2alias="sim",
    polarize=TRUE, 
    outgroupSeqs = list(risa_popIDs[["popDat_justSim"]][["yak"]]))


## popDat_bothSpecies
risa_MKresults[["popDat_bothSpecies_unpolarizedA"]] <- doMKtest(
    risaFiles[["popDat_bothSpecies"]], 
    pop1seqs=risa_popIDs[["popDat_bothSpecies"]][["sim"]], pop1alias="sim",
    pop2seqs=risa_popIDs[["popDat_bothSpecies"]][["mel"]], pop2alias="mel")

risa_MKresults[["popDat_bothSpecies_unpolarizedB"]] <- doMKtest(
    risaFiles[["popDat_bothSpecies"]], 
    pop1seqs=risa_popIDs[["popDat_bothSpecies"]][["mel"]], pop1alias="mel",
    pop2seqs=risa_popIDs[["popDat_bothSpecies"]][["sim"]], pop2alias="sim")

risa_MKresults[["popDat_bothSpecies_polarizedA"]] <- doMKtest(
    risaFiles[["popDat_bothSpecies"]], 
    pop1seqs=risa_popIDs[["popDat_bothSpecies"]][["sim"]], pop1alias="sim",
    pop2seqs=risa_popIDs[["popDat_bothSpecies"]][["mel"]], pop2alias="mel",
    polarize=TRUE, 
    outgroupSeqs = list(risa_popIDs[["popDat_bothSpecies"]][["yak"]]))

risa_MKresults[["popDat_bothSpecies_polarizedB"]] <- doMKtest(
    risaFiles[["popDat_bothSpecies"]], 
    pop1seqs=risa_popIDs[["popDat_bothSpecies"]][["mel"]], pop1alias="mel",
    pop2seqs=risa_popIDs[["popDat_bothSpecies"]][["sim"]], pop2alias="sim",
    polarize=TRUE, 
    outgroupSeqs = list(risa_popIDs[["popDat_bothSpecies"]][["yak"]]))


risa_MKresults[["popDat_bothSpec_polB_noRare"]] <- doMKtest(
    risaFiles[["popDat_bothSpecies"]], 
    pop1seqs=risa_popIDs[["popDat_bothSpecies"]][["mel"]], pop1alias="mel",
    pop2seqs=risa_popIDs[["popDat_bothSpecies"]][["sim"]], pop2alias="sim",
    polarize=TRUE, 
    outgroupSeqs = list(risa_popIDs[["popDat_bothSpecies"]][["yak"]]),
    filterRareAlleles = TRUE, alleleFreqThreshold=0.05)



## combine and save to Excel files
# just unpolarized
risa_MKresults_unpolarized <- combineMKresults(risa_MKresults[ grep("_unpolarized", names(risa_MKresults)) ], 
                                               outFile="risa_MKresults_unpolarized.xlsx",
                                               outDir=here("data/Risa/2022_Dec7"),
                                               getGeneNames=FALSE)
# just polarized
risa_MKresults_polarized <- combineMKresults(risa_MKresults[ grep("_polarized", names(risa_MKresults)) ], 
                                             outFile="risa_MKresults_polarized.xlsx",
                                             outDir=here("data/Risa/2022_Dec7"),
                                             getGeneNames=FALSE, 
                                             keepNA=FALSE, NAcharacter="")

# all of them together
risa_MKresults_all <- combineMKresults(risa_MKresults, 
                                       outFile="risa_MKresults_all_addColumns.xlsx",
                                       outDir=here("data/Risa/2022_Dec7"),
                                       getGeneNames=FALSE)

### sanity check:  now, after fixing code:
# - results are symmetrical, no matter which population is pop1 and which is pop2
# - I do think results are making more sense

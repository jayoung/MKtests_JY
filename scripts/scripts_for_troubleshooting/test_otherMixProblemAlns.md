test_otherMiscProblemAlns.Rmd
================
Janet Young

2025-02-07

# deal with ragged alignments

i.e. sequences are not all the same length

I have now built in an error check for this

``` r
raggedAlnFile <- here("test_data/testData_for_troubleshooting/test_otherMiscProblemAlns/test_raggedAln.fa")

## get IDs of individuals in each population
raggedAln <- readBStringSet(raggedAlnFile)
names(raggedAln) <- sapply(strsplit(names(raggedAln), " "), "[[", 1)
raggedAlnSpecies <- gsub("\\d","",names(raggedAln))
raggedAlnPopulations <- split(names(raggedAln), raggedAlnSpecies)

# test
MKresults <- doMKtest(raggedAlnFile, 
                      pop1seqs=raggedAlnPopulations[["pongo"]],
                      pop2seqs=raggedAlnPopulations[["trachy"]])
```

    ## 
    ## ##### reading alignment from file /fh/fast/malik_h/user/jayoung/MKtest/MKtests_JY/test_data/testData_for_troubleshooting/test_otherMiscProblemAlns/test_raggedAln.fa

    ## Error in doMKtest(raggedAlnFile, pop1seqs = raggedAlnPopulations[["pongo"]], : 
    ## 
    ## ERROR - your alignment looks bad. Sequences are not all the same length as each other

# alignments that are not an even multiple of 3 in length

I check for that.

``` r
partialCodonAlnFile <- here("test_data/testData_for_troubleshooting/test_otherMiscProblemAlns/test_partialCodonAln.fa")

MKresults <- doMKtest(partialCodonAlnFile, 
                      pop1seqs=raggedAlnPopulations[["pongo"]],
                      pop2seqs=raggedAlnPopulations[["trachy"]])
```

    ## 
    ## ##### reading alignment from file /fh/fast/malik_h/user/jayoung/MKtest/MKtests_JY/test_data/testData_for_troubleshooting/test_otherMiscProblemAlns/test_partialCodonAln.fa

    ## Error in doMKtest(partialCodonAlnFile, pop1seqs = raggedAlnPopulations[["pongo"]], : ERROR - alignment length is not a multiple of 3: 4

# Check whether it works with mixed upper and lower case seqs

``` r
upperLowerAlnFile <- here("test_data/testData_for_troubleshooting/test_otherMiscProblemAlns/MKTwebsite_testAln_mixUpperLowerCase.fa")

# get IDs of individuals in each population
upperLowerAln <- readBStringSet(upperLowerAlnFile)
names(upperLowerAln) <- sapply(strsplit(names(upperLowerAln), " "), "[[", 1)
upperLowerAlnSpecies <- gsub("\\d","",names(upperLowerAln))
upperLowerAlnPopulations <- split(names(upperLowerAln), upperLowerAlnSpecies)

MKresults <- doMKtest(upperLowerAlnFile, 
                      pop1seqs=upperLowerAlnPopulations[["pongo"]],
                      pop2seqs=upperLowerAlnPopulations[["trachy"]])
```

    ## 
    ## ##### reading alignment from file /fh/fast/malik_h/user/jayoung/MKtest/MKtests_JY/test_data/testData_for_troubleshooting/test_otherMiscProblemAlns/MKTwebsite_testAln_mixUpperLowerCase.fa 
    ##     splitting alignment into groups and tabulating nucleotides

    ## Warning in doMKtest(upperLowerAlnFile, pop1seqs = upperLowerAlnPopulations[["pongo"]], : this alignment has a codon at the end with only stop codons - we will strip it out and not count any changes in this codon

    ##     looking at allele frequencies
    ##     inferring ancestors
    ##     starting output table
    ##     categorizing population 1 vs 2 fixed changes
    ##     categorizing polymorphisms
    ##     doing MK tests
    ##         seqs_not_used

# check whether population name lists are overlapping

user might accidentally specify populations with names in common.
that’ll mess things up.

``` r
MKresults <- doMKtest(upperLowerAlnFile, 
                      pop1seqs=upperLowerAlnPopulations[["pongo"]],
                      pop2seqs=upperLowerAlnPopulations[["pongo"]])
```

    ## 
    ## ##### reading alignment from file /fh/fast/malik_h/user/jayoung/MKtest/MKtests_JY/test_data/testData_for_troubleshooting/test_otherMiscProblemAlns/MKTwebsite_testAln_mixUpperLowerCase.fa

    ## Error in doMKtest(upperLowerAlnFile, pop1seqs = upperLowerAlnPopulations[["pongo"]], : 
    ## 
    ## ERROR - there are sequence names found in both pop1seqs and pop2seqs - that's not right.
    ## The overlapping names are: pongo1 pongo2

# check for empty alignment file

``` r
emptyAlnFile <- here("test_data/testData_for_troubleshooting/test_otherMiscProblemAlns/test_emptyAlnFile.fa")

MKresults <- doMKtest(emptyAlnFile, 
                      pop1seqs=raggedAlnPopulations[["pongo"]],
                      pop2seqs=raggedAlnPopulations[["trachy"]])
```

    ## 
    ## ##### reading alignment from file /fh/fast/malik_h/user/jayoung/MKtest/MKtests_JY/test_data/testData_for_troubleshooting/test_otherMiscProblemAlns/test_emptyAlnFile.fa

    ## Error in doMKtest(emptyAlnFile, pop1seqs = raggedAlnPopulations[["pongo"]], : 
    ## 
    ## ERROR - the alignment file you specified exists but is empty: /fh/fast/malik_h/user/jayoung/MKtest/MKtests_JY/test_data/testData_for_troubleshooting/test_otherMiscProblemAlns/test_emptyAlnFile.fa

# check for stop codons fixed in one population but not the other (or in the outgroup)

``` r
stopCodonAllSeqs <- here("test_data/testData_for_troubleshooting/test_otherMiscProblemAlns/MKTwebsite_testAln.fa")
stopCodonAln1 <- here("test_data/testData_for_troubleshooting/test_otherMiscProblemAlns/test_onePopHasStop.fa")
stopCodonAln2 <- here("test_data/testData_for_troubleshooting/test_otherMiscProblemAlns/test_outgroupHasStop.fa")

# test
MKresults <- doMKtest(stopCodonAllSeqs, 
                      pop1seqs=raggedAlnPopulations[["pongo"]],
                      pop2seqs=raggedAlnPopulations[["trachy"]])
```

    ## 
    ## ##### reading alignment from file /fh/fast/malik_h/user/jayoung/MKtest/MKtests_JY/test_data/testData_for_troubleshooting/test_otherMiscProblemAlns/MKTwebsite_testAln.fa 
    ##     splitting alignment into groups and tabulating nucleotides

    ## Warning in doMKtest(stopCodonAllSeqs, pop1seqs = raggedAlnPopulations[["pongo"]], : this alignment has a codon at the end with only stop codons - we will strip it out and not count any changes in this codon

    ##     looking at allele frequencies
    ##     inferring ancestors
    ##     starting output table
    ##     categorizing population 1 vs 2 fixed changes
    ##     categorizing polymorphisms
    ##     doing MK tests
    ##         seqs_not_used

``` r
MKresults <- doMKtest(stopCodonAln1, 
                      pop1seqs=raggedAlnPopulations[["pongo"]],
                      pop2seqs=raggedAlnPopulations[["trachy"]])
```

    ## 
    ## ##### reading alignment from file /fh/fast/malik_h/user/jayoung/MKtest/MKtests_JY/test_data/testData_for_troubleshooting/test_otherMiscProblemAlns/test_onePopHasStop.fa 
    ##     splitting alignment into groups and tabulating nucleotides
    ##     looking at allele frequencies
    ##     inferring ancestors
    ## all codons TGA 
    ## non-stop codons

    ## Error in FUN(X[[i]], ...): 
    ## 
    ## ERROR in population pop2_anc, in codon 133 - all codon possibilities are stop codons - that's odd

# work on combining outputs

I want to be able to combine polarized and unpolarized results into the
same output file

``` r
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
```

    ## 
    ## ##### reading alignment from file /fh/fast/malik_h/user/jayoung/MKtest/MKtests_JY/data/Risa/2022_Dec7/abo_simpop_mel_yak.fasta 
    ##     splitting alignment into groups and tabulating nucleotides
    ##     looking at allele frequencies
    ##     inferring ancestors
    ##     starting output table
    ##     categorizing population 1 vs 2 fixed changes
    ##     categorizing polymorphisms
    ##     doing MK tests
    ##         seqs_not_used D.yakuba

``` r
risa_MKresults[["popDat_justSim_unpolarizedB"]] <- doMKtest(
    risaFiles[["popDat_justSim"]], 
    pop1seqs=risa_popIDs[["popDat_justSim"]][["mel"]], pop1alias="mel",
    pop2seqs=risa_popIDs[["popDat_justSim"]][["sim"]], pop2alias="sim")
```

    ## 
    ## ##### reading alignment from file /fh/fast/malik_h/user/jayoung/MKtest/MKtests_JY/data/Risa/2022_Dec7/abo_simpop_mel_yak.fasta 
    ##     splitting alignment into groups and tabulating nucleotides
    ##     looking at allele frequencies
    ##     inferring ancestors
    ##     starting output table
    ##     categorizing population 1 vs 2 fixed changes
    ##     categorizing polymorphisms
    ##     doing MK tests
    ##         seqs_not_used D.yakuba

``` r
risa_MKresults[["popDat_justSim_polarizedA"]] <- doMKtest(
    risaFiles[["popDat_justSim"]], 
    pop1seqs=risa_popIDs[["popDat_justSim"]][["sim"]], pop1alias="sim",
    pop2seqs=risa_popIDs[["popDat_justSim"]][["mel"]], pop2alias="mel",
    polarize=TRUE, 
    outgroupSeqs = list(risa_popIDs[["popDat_justSim"]][["yak"]]))
```

    ## 
    ## ##### reading alignment from file /fh/fast/malik_h/user/jayoung/MKtest/MKtests_JY/data/Risa/2022_Dec7/abo_simpop_mel_yak.fasta 
    ##     splitting alignment into groups and tabulating nucleotides
    ##     looking at allele frequencies
    ##     inferring ancestors
    ##     starting output table
    ##     categorizing population 1 vs 2 fixed changes
    ##     categorizing polymorphisms
    ##     getting polarized changes
    ## 
    ##     doing MK tests
    ##         seqs_not_used

``` r
risa_MKresults[["popDat_justSim_polarizedB"]] <- doMKtest(
    risaFiles[["popDat_justSim"]], 
    pop1seqs=risa_popIDs[["popDat_justSim"]][["mel"]], pop1alias="mel",
    pop2seqs=risa_popIDs[["popDat_justSim"]][["sim"]], pop2alias="sim",
    polarize=TRUE, 
    outgroupSeqs = list(risa_popIDs[["popDat_justSim"]][["yak"]]))
```

    ## 
    ## ##### reading alignment from file /fh/fast/malik_h/user/jayoung/MKtest/MKtests_JY/data/Risa/2022_Dec7/abo_simpop_mel_yak.fasta 
    ##     splitting alignment into groups and tabulating nucleotides
    ##     looking at allele frequencies
    ##     inferring ancestors
    ##     starting output table
    ##     categorizing population 1 vs 2 fixed changes
    ##     categorizing polymorphisms
    ##     getting polarized changes
    ## 
    ##     doing MK tests
    ##         seqs_not_used

``` r
## popDat_bothSpecies
risa_MKresults[["popDat_bothSpecies_unpolarizedA"]] <- doMKtest(
    risaFiles[["popDat_bothSpecies"]], 
    pop1seqs=risa_popIDs[["popDat_bothSpecies"]][["sim"]], pop1alias="sim",
    pop2seqs=risa_popIDs[["popDat_bothSpecies"]][["mel"]], pop2alias="mel")
```

    ## 
    ## ##### reading alignment from file /fh/fast/malik_h/user/jayoung/MKtest/MKtests_JY/data/Risa/2022_Dec7/abo_melpop_simpop_yak.fasta 
    ##     splitting alignment into groups and tabulating nucleotides
    ##     looking at allele frequencies
    ##     inferring ancestors
    ##     starting output table
    ##     categorizing population 1 vs 2 fixed changes
    ##     categorizing polymorphisms
    ##     doing MK tests
    ##         seqs_not_used D.yakuba

``` r
risa_MKresults[["popDat_bothSpecies_unpolarizedB"]] <- doMKtest(
    risaFiles[["popDat_bothSpecies"]], 
    pop1seqs=risa_popIDs[["popDat_bothSpecies"]][["mel"]], pop1alias="mel",
    pop2seqs=risa_popIDs[["popDat_bothSpecies"]][["sim"]], pop2alias="sim")
```

    ## 
    ## ##### reading alignment from file /fh/fast/malik_h/user/jayoung/MKtest/MKtests_JY/data/Risa/2022_Dec7/abo_melpop_simpop_yak.fasta 
    ##     splitting alignment into groups and tabulating nucleotides
    ##     looking at allele frequencies
    ##     inferring ancestors
    ##     starting output table
    ##     categorizing population 1 vs 2 fixed changes
    ##     categorizing polymorphisms
    ##     doing MK tests
    ##         seqs_not_used D.yakuba

``` r
risa_MKresults[["popDat_bothSpecies_polarizedA"]] <- doMKtest(
    risaFiles[["popDat_bothSpecies"]], 
    pop1seqs=risa_popIDs[["popDat_bothSpecies"]][["sim"]], pop1alias="sim",
    pop2seqs=risa_popIDs[["popDat_bothSpecies"]][["mel"]], pop2alias="mel",
    polarize=TRUE, 
    outgroupSeqs = list(risa_popIDs[["popDat_bothSpecies"]][["yak"]]))
```

    ## 
    ## ##### reading alignment from file /fh/fast/malik_h/user/jayoung/MKtest/MKtests_JY/data/Risa/2022_Dec7/abo_melpop_simpop_yak.fasta 
    ##     splitting alignment into groups and tabulating nucleotides
    ##     looking at allele frequencies
    ##     inferring ancestors
    ##     starting output table
    ##     categorizing population 1 vs 2 fixed changes
    ##     categorizing polymorphisms
    ##     getting polarized changes
    ## 
    ##     doing MK tests
    ##         seqs_not_used

``` r
risa_MKresults[["popDat_bothSpecies_polarizedB"]] <- doMKtest(
    risaFiles[["popDat_bothSpecies"]], 
    pop1seqs=risa_popIDs[["popDat_bothSpecies"]][["mel"]], pop1alias="mel",
    pop2seqs=risa_popIDs[["popDat_bothSpecies"]][["sim"]], pop2alias="sim",
    polarize=TRUE, 
    outgroupSeqs = list(risa_popIDs[["popDat_bothSpecies"]][["yak"]]))
```

    ## 
    ## ##### reading alignment from file /fh/fast/malik_h/user/jayoung/MKtest/MKtests_JY/data/Risa/2022_Dec7/abo_melpop_simpop_yak.fasta 
    ##     splitting alignment into groups and tabulating nucleotides
    ##     looking at allele frequencies
    ##     inferring ancestors
    ##     starting output table
    ##     categorizing population 1 vs 2 fixed changes
    ##     categorizing polymorphisms
    ##     getting polarized changes
    ## 
    ##     doing MK tests
    ##         seqs_not_used

``` r
## combine and save to Excel files

## without excel file
risa_MKresults_unpolarized <- combineMKresults(risa_MKresults[ grep("_unpolarized", names(risa_MKresults)) ])
risa_MKresults_polarized <- combineMKresults(risa_MKresults[ grep("_polarized", names(risa_MKresults)) ])


### with excel file
risa_MKresults_unpolarized <- combineMKresults(
    risa_MKresults[ grep("_unpolarized", names(risa_MKresults)) ],
    outFile="risa_MKresults_unpolarized_v1.xlsx",
    outDir=here("test_data/testData_for_troubleshooting/test_otherMiscProblemAlns"))
```

    ## Loading required package: openxlsx

``` r
### with excel file keepNA=FALSE, NAcharacter="", 
risa_MKresults_unpolarized <- combineMKresults(
    risa_MKresults[ grep("_unpolarized", names(risa_MKresults)) ],
    outFile="risa_MKresults_unpolarized_v2.xlsx",
    outDir=here("test_data/testData_for_troubleshooting/test_otherMiscProblemAlns"),
    keepNA=FALSE, NAcharacter="")

risa_MKresults_all <- combineMKresults(risa_MKresults, extraVerbose=TRUE)
```

    ## colnames are now input num_seqs num_seqs_pop1 num_seqs_pop2 num_seqs_outgroup seqs_not_used length_NT length_AA region_start_coord_codons region_end_coord_codons filter_rare_alleles rare_allele_freq_threshold parameter_combining_approach outgroup first_pop1_seq first_pop2_seq pop1_vs_pop2_chiSq_pVal pop1_vs_pop2_FET_pVal pop1_vs_pop2_Dn pop1_vs_pop2_Ds pop1_vs_pop2_Pn pop1_vs_pop2_Ps pop1_vs_pop2_NI pop1_vs_pop2_alpha pop1_vs_pop2_result polarized_pop1_chiSq_pVal polarized_pop1_FET_pVal polarized_pop1_Dn polarized_pop1_Ds pop1_Pn pop1_Ps polarized_pop1_NI polarized_pop1_alpha polarized_pop1_result polarized_pop2_chiSq_pVal polarized_pop2_FET_pVal polarized_pop2_Dn polarized_pop2_Ds pop2_Pn pop2_Ps polarized_pop2_NI polarized_pop2_alpha polarized_pop2_result total_polarized_Dn total_polarized_Ds

``` r
risa_MKresults_all <- combineMKresults(
    risa_MKresults, 
    outFile="risa_MKresults_all.xlsx",
    outDir=here("test_data/testData_for_troubleshooting/test_otherMiscProblemAlns"), 
    extraVerbose=TRUE)
```

    ## colnames are now input num_seqs num_seqs_pop1 num_seqs_pop2 num_seqs_outgroup seqs_not_used length_NT length_AA region_start_coord_codons region_end_coord_codons filter_rare_alleles rare_allele_freq_threshold parameter_combining_approach outgroup first_pop1_seq first_pop2_seq pop1_vs_pop2_chiSq_pVal pop1_vs_pop2_FET_pVal pop1_vs_pop2_Dn pop1_vs_pop2_Ds pop1_vs_pop2_Pn pop1_vs_pop2_Ps pop1_vs_pop2_NI pop1_vs_pop2_alpha pop1_vs_pop2_result polarized_pop1_chiSq_pVal polarized_pop1_FET_pVal polarized_pop1_Dn polarized_pop1_Ds pop1_Pn pop1_Ps polarized_pop1_NI polarized_pop1_alpha polarized_pop1_result polarized_pop2_chiSq_pVal polarized_pop2_FET_pVal polarized_pop2_Dn polarized_pop2_Ds pop2_Pn pop2_Ps polarized_pop2_NI polarized_pop2_alpha polarized_pop2_result total_polarized_Dn total_polarized_Ds 
    ## here1 colnames are now gene_name input num_seqs num_seqs_pop1 num_seqs_pop2 num_seqs_outgroup seqs_not_used length_NT length_AA region_start_coord_codons region_end_coord_codons filter_rare_alleles rare_allele_freq_threshold parameter_combining_approach outgroup first_pop1_seq first_pop2_seq pop1_vs_pop2_chiSq_pVal pop1_vs_pop2_FET_pVal pop1_vs_pop2_Dn pop1_vs_pop2_Ds pop1_vs_pop2_Pn pop1_vs_pop2_Ps pop1_vs_pop2_NI pop1_vs_pop2_alpha pop1_vs_pop2_result polarized_pop1_chiSq_pVal polarized_pop1_FET_pVal polarized_pop1_Dn polarized_pop1_Ds pop1_Pn pop1_Ps polarized_pop1_NI polarized_pop1_alpha polarized_pop1_result polarized_pop2_chiSq_pVal polarized_pop2_FET_pVal polarized_pop2_Dn polarized_pop2_Ds pop2_Pn pop2_Ps polarized_pop2_NI polarized_pop2_alpha polarized_pop2_result total_polarized_Dn total_polarized_Ds 
    ## here2 colnames are now gene name input num seqs num seqs pop1 num seqs pop2 num seqs outgroup seqs not used length NT length AA region start coord codons region end coord codons filter rare alleles rare allele freq threshold parameter combining approach outgroup first pop1 seq first pop2 seq pop1 vs pop2 chiSq pVal pop1 vs pop2 FET pVal pop1 vs pop2 Dn pop1 vs pop2 Ds pop1 vs pop2 Pn pop1 vs pop2 Ps pop1 vs pop2 NI pop1 vs pop2 alpha pop1 vs pop2 result polarized pop1 chiSq pVal polarized pop1 FET pVal polarized pop1 Dn polarized pop1 Ds pop1 Pn pop1 Ps polarized pop1 NI polarized pop1 alpha polarized pop1 result polarized pop2 chiSq pVal polarized pop2 FET pVal polarized pop2 Dn polarized pop2 Ds pop2 Pn pop2 Ps polarized pop2 NI polarized pop2 alpha polarized pop2 result total polarized Dn total polarized Ds 
    ## columnsToRound  18 19 24 25 27 28 33 34 36 37 42 43

``` r
### weird - ALL OF the column headers are blank with this:
risa_MKresults_all_2 <- combineMKresults(
    risa_MKresults, 
    outFile="risa_MKresults_all_v2.xlsx",
    outDir=here("test_data/testData_for_troubleshooting/test_otherMiscProblemAlns"),
    keepNA=FALSE, NAcharacter="", 
    extraVerbose=TRUE)
```

    ## Error in FUN(X[[i]], ...): 
    ## 
    ## ERROR - haven't figured out the correct way to code this situation. We are trying to merge unpolarized and polarized results into a single output file, and we are trying to use the keepNA=FALSE option while exporting to Excel. It's acting weirdly and I haven't done figured out why yet

``` r
### weird - ALL OF the column headers are are blank using this:
risa_MKresults_all_3 <- combineMKresults(
    risa_MKresults, 
    outFile="risa_MKresults_all_v3.xlsx",
    outDir=here("test_data/testData_for_troubleshooting/test_otherMiscProblemAlns"),
    keepNA=FALSE, NAcharacter="N.A.", 
    extraVerbose=TRUE)
```

    ## Error in FUN(X[[i]], ...): 
    ## 
    ## ERROR - haven't figured out the correct way to code this situation. We are trying to merge unpolarized and polarized results into a single output file, and we are trying to use the keepNA=FALSE option while exporting to Excel. It's acting weirdly and I haven't done figured out why yet

# Finished

``` r
sessionInfo()
```

    ## R version 4.4.0 (2024-04-24)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 18.04.6 LTS
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: FlexiBLAS OPENBLAS;  LAPACK version 3.11.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: America/Los_Angeles
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] openxlsx_4.2.5.2    Biostrings_2.72.0   GenomeInfoDb_1.40.0
    ##  [4] XVector_0.44.0      IRanges_2.38.0      S4Vectors_0.42.0   
    ##  [7] BiocGenerics_0.50.0 kableExtra_1.4.0    lubridate_1.9.3    
    ## [10] forcats_1.0.0       stringr_1.5.1       dplyr_1.1.4        
    ## [13] purrr_1.0.2         readr_2.1.5         tidyr_1.3.1        
    ## [16] tibble_3.2.1        ggplot2_3.5.1       tidyverse_2.0.0    
    ## [19] here_1.0.1         
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] utf8_1.2.4              generics_0.1.3          xml2_1.3.6             
    ##  [4] stringi_1.8.4           hms_1.1.3               digest_0.6.35          
    ##  [7] magrittr_2.0.3          evaluate_0.23           grid_4.4.0             
    ## [10] timechange_0.3.0        fastmap_1.2.0           jsonlite_1.8.8         
    ## [13] rprojroot_2.0.4         zip_2.3.1               httr_1.4.7             
    ## [16] fansi_1.0.6             UCSC.utils_1.0.0        viridisLite_0.4.2      
    ## [19] scales_1.3.0            cli_3.6.2               crayon_1.5.2           
    ## [22] rlang_1.1.4             munsell_0.5.1           withr_3.0.0            
    ## [25] yaml_2.3.8              tools_4.4.0             tzdb_0.4.0             
    ## [28] colorspace_2.1-0        GenomeInfoDbData_1.2.12 vctrs_0.6.5            
    ## [31] R6_2.5.1                lifecycle_1.0.4         zlibbioc_1.50.0        
    ## [34] pkgconfig_2.0.3         pillar_1.9.0            gtable_0.3.5           
    ## [37] Rcpp_1.0.12             glue_1.7.0              systemfonts_1.1.0      
    ## [40] xfun_0.44               tidyselect_1.2.1        rstudioapi_0.16.0      
    ## [43] knitr_1.46              htmltools_0.5.8.1       rmarkdown_2.26         
    ## [46] svglite_2.1.3           compiler_4.4.0

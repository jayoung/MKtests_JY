### https://en.wikipedia.org/wiki/McDonald%E2%80%93Kreitman_test

require(Biostrings)
source("scripts/MKfunctions_utilities.R")
source("scripts/MKfunctions_plottingAndExcelOutput.R")
load("scripts/codonPaths/codonPathsFromBioperl.Rdata")

## on fast the base dir is /fh/fast/malik_h/user/jayoung/MKtest/MKtests_JY
## on mac (local) the base dir is /Users/jayoung/Desktop/mac_workStuff/mac_MKtests/MKtests_JY
## on mac (server) the base dir is /Volumes/malik_h/user/jayoung/MKtest/MKtests_JY



#### reconstructAncestor: a function to reconstruct ancestors using nucleotide count tables. Useful for polarized tests, but perhaps also just to get mel ancestor using sim as outgroup.
## arguments:
# ingroupName(s) e.g. "mel" or "sim"
# orderedOutgroupNames e.g. "sim" or c("sim", "out")

reconstructAncestor <- function(aln_tables, aln_df, fixed_nucs, seqGroups, ingroupNames, orderedOutgroupNames) {
    numNT <- dim(aln_tables[[1]])[2]
    anc <- rep("?", numNT)
    
    if (sum(!ingroupNames %in% names(aln_tables))>0) {
        stop("ERROR - one or more unrecognized ingroup name:",ingroupNames,"\n\n")
    }
    if (sum(!orderedOutgroupNames %in% names(aln_tables))>0) {
        stop("ERROR - one or more unrecognized outgroup name:", orderedOutgroupNames, "\n\n")
    }
    # if there is more than one ingroup I want to add their count tables and seqnames together. I also want to rederive the vector that tells me whether a position is fixed
    ingroup_table <- aln_tables[[ ingroupNames[1] ]]
    fixedPositions <- which(fixed_nucs[[ ingroupNames[1] ]])
    
    if (length(ingroupNames)>1) {
        for (thisIngroup in ingroupNames[2:length(ingroupNames)]) {
            ingroup_table <- ingroup_table + aln_tables[[ thisIngroup ]]
        }
        fixedPositions <- apply(ingroup_table, 2, function(y) {
            sum(y>0) == 1
        })
        names(fixedPositions) <- NULL
        fixedPositions <- which(fixedPositions)
    }
    
    # for positions that are fixed, I take that nucleotide in the first seq of the ingroup
    firstIngroupSeq <- seqGroups[[ ingroupNames[1] ]][1]
    anc[ fixedPositions ] <- as.character(aln_df[ firstIngroupSeq,fixedPositions])
    
    # for positions that are polymorphic, I use the outgroup(s)
    # use each outgroup in turn to assign ancestors at any remaining uncertain positions
    for (thisOutgroup in orderedOutgroupNames) {
        #cat("        Using",thisOutgroup,"seq(s) to reconstruct ancestor of",ingroupNames,"\n")
        #cat("            there are",sum(anc == "?"),"uncertain positions\n") 
        if (sum(anc == "?")>0) {
            polymorphicPositions <- which(anc == "?")
            tempAncs <- lapply(polymorphicPositions, function(i) {
                ingroupCounts <- ingroup_table[,i]
                ingroupNucs <- rownames( ingroup_table ) [ which(ingroupCounts>0) ]
                outgroupCounts <- aln_tables[[thisOutgroup]][,i]
                outgroupNucs <- rownames( aln_tables[[thisOutgroup]] ) [ which(outgroupCounts>0) ]
                # if only 1 of the ingroup nucs is in the outgroup nucs, that must be the ancestor.  Other possibilities: 0 or >1
                possibleAncs <- intersect(ingroupNucs, outgroupNucs)
                if(length(possibleAncs)==1) {
                    return(possibleAncs[1])
                } else {
                    return("?")
                }
            })
            anc[polymorphicPositions] <- unlist(tempAncs)
        }
        #cat("            there are now",sum(anc == "?"),"uncertain positions\n") 
    }
    return(anc)
}


####### reconstructAncestor_includeUncertainty: a function to reconstruct ancestors at individual positions (NOT codons), now including uncertainty
## uncertainty is only regarding ambiguous reconstructions. If any of the extant sequences have N we ignore the N
# ingroup is a list like the output of alnSlicesUniqueSeqs for the ingroup we are trying to get ancestor for
# outgroup(s) is a list of lists, for outgroups in order. We use the first outgroup to resolve uncertainty if possible, if not we try the second, etc
reconstructAncestor_includeUncertainty <- function(ingroup, outgroups,
                                                   extraVerbose=FALSE) {
    if (nchar(ingroup[[1]][1]) != 1) {
        stop("ERROR - alignment should be split into single nucleotide slices for the reconstructAncestor_includeUncertainty function\n\n")
    }
    numPositions <- length(ingroup)
    ancestor <- lapply(1:numPositions, function(i) {
        ingroupThisPos <- ingroup[[i]]
        if (length(ingroupThisPos)==1) {return(ingroupThisPos)}
        outgroupsThisPos <- lapply(outgroups, "[[", i)
        if(extraVerbose) {
            cat("\n            ingroupThisPos:\n")
            print(ingroupThisPos)
            cat("\n            outgroupsThisPos:\n")
            print(outgroupsThisPos)
            cat("\n")
        }
        # try each outgroup in turn for a unique solution
        for (outgroupThisPos in outgroupsThisPos) {
            sharedNucs <- intersect(ingroupThisPos,outgroupThisPos)
            if (length(sharedNucs)==1) {
                return(sharedNucs) # goes to the next position
            }
        }
        # if we still have not worked it out:
        return(ingroupThisPos)
    })
    if(extraVerbose) {
        cat("            ancestor:\n")
        print(ancestor)
        cat("\n")
    }
    return(ancestor)
}

#### get counts of synon and non-synon changes between sets of codons (character vectors)
getCodonChangeCounts <- function(codonsBefore, codonsAfter, 
                                 paths=codonPaths, 
                                 extraVerbose=FALSE) {

    if (length(codonsBefore) != length(codonsAfter)){
        stop("ERROR - codonsBefore and codonsAfter are different lengths\n\n")    
    }
    if (extraVerbose) {
        cat("codonsBefore", codonsBefore,"\n")
        cat("codonsAfter", codonsAfter,"\n")
    }
    ## make sure we have at least one non-stop codon
    outputIfWeFindStopCodons <- data.frame(codon=1:length(codonsBefore), 
                                           before=codonsBefore, 
                                           after=codonsAfter,
                                           ns=0,
                                           s=0,
                                           tot=0 )
    codonsBeforeWithoutStops <- setdiff(codonsBefore, stopCodons)
    codonsAfterWithoutStops <- setdiff(codonsAfter, stopCodons)
    if(length(codonsBeforeWithoutStops)==0 | length(codonsAfterWithoutStops)==0 ) {
        cat("\n        NOTE from getCodonChangeCounts function. There was a weird situation - a combination of ancestral allele+polymorphic allele within a codon that would create a stop codon. Unlikely to exist in nature. This situation tends to arise when a polymorphism arises within a codon where there is >1 change from the ancestral sequence. OK to ignore this warning\n\n")
        return(outputIfWeFindStopCodons)
    }
    
    nsCountsEachCodon <- rep(NA,length(codonsBefore))
    sCountsEachCodon <- rep(NA,length(codonsBefore))
    
    ## the easy case is identical codons (even if they contain - but if one contains "?" I will NOT call them identical)
    identicalCodons <- which(codonsBefore==codonsAfter & 
                                 !grepl("?", codonsBefore, fixed=TRUE)  & 
                                 !grepl("?", codonsAfter, fixed=TRUE) )
    nsCountsEachCodon[identicalCodons] <- 0
    sCountsEachCodon[identicalCodons] <- 0
    codonsJoined <- paste(codonsBefore, codonsAfter, sep="")
    
    ## codons that differ and are unambiguous, and we can get counts for those
    differentCodons <- which(codonsBefore!=codonsAfter & 
                                 !grepl("[-N\\?]", codonsJoined, perl=TRUE, ignore.case=TRUE))
    if(extraVerbose) {
        cat("checking codonsJoined ",codonsJoined,"\n\n")
    }
    # make sure all codons are in the table
    checkCodons <- unique(codonsJoined[differentCodons])
    if (sum(!checkCodons %in% paths[,"path"])>0) {
        problemCodonPositions <- which(!checkCodons %in% paths[,"path"])
        missingCodons <- checkCodons[which(!checkCodons %in% paths[,"path"])]
        stop("ERROR - found a codon combination missing from the paths table, probably contains a stop codon: ",missingCodons," at position ",problemCodonPositions,"\n\n")
    }
    
    pathsRows <- match( codonsJoined[differentCodons], paths[,"path"] )
    nsCountsEachCodon[differentCodons] <- paths[pathsRows,"ns"]
    sCountsEachCodon[differentCodons] <- paths[pathsRows,"s"]
    output <- data.frame(codon=1:length(codonsBefore), 
                         before=codonsBefore, 
                         after=codonsAfter,
                         ns=nsCountsEachCodon,
                         s=sCountsEachCodon,
                         tot= (nsCountsEachCodon+sCountsEachCodon) )
    
    return(output)
    
    ## xx maybe do something with codons that contain ?
    # I could perhaps average out over the uncertainty ?
}

#### getCodonChangeCountsFromCodonList: a function to get counts of synon and non-synon changes between sets of codons (character vectors)
# this version allows for ambiguity within a codon, and averages result over that ambiguity
# now, input= two lists to compare. Each list has one element per codon, containing character vector of the codons to compare
# combiningApproach can be:
#     "mean" (take mean ns and mean s over all possible combinations) or 
#     "conservative" (take the count combination that minimizes ns)
getCodonChangeCountsFromCodonList <- function(codonListBefore, codonListAfter, 
                                              paths=codonPaths, 
                                              combiningApproach="conservative",
                                              extraVerbose=FALSE) {
    ### some checks
    if (length(codonListBefore) != length(codonListAfter)){
        stop("\n\nERROR in getCodonChangeCountsFromCodonList - codonListBefore and codonListAfter are different lengths\n\n")    
    }
    if (!combiningApproach %in% c("mean","conservative","conservativeOld") ) {
        stop("\n\nERROR in getCodonChangeCountsFromCodonList - the combiningApproach option can only be mean or conservative or conservativeOld\n\n")
    }
    
    if(sum(stopCodons %in% codonListBefore) > 0) {
        stop("\n\nERRR in getCodonChangeCountsFromCodonList - there are stop codon(s) in codonListBefore\n\n")
    }
    if(sum(stopCodons %in% codonListAfter) > 0) {
        stop("\n\nERRR in getCodonChangeCountsFromCodonList - there are stop codon(s) in codonListAfter\n\n")
    }
    
    ### get counts, one codon position at a time
    numCodons <- length(codonListBefore)
    counts <- lapply(1:numCodons, function(i) {
        before <- codonListBefore[[i]]
        after <- codonListAfter[[i]]
        # generate all possible comparisons between before and after
        numComparisons <- length(before) * length(after)
        beforeAll <- rep( before, length(after))
        afterAll <- rep( after, each=length(before))
        # each possible combination shows up as a row in this result table
        results <- getCodonChangeCounts(beforeAll,afterAll, extraVerbose=extraVerbose)
        if(extraVerbose) {
            cat("            results BEFORE combining:\n")
            print(results)
            cat("\n")
        }
        # now get a single ns and s considering all possible combinations
        if (combiningApproach=="mean") {
            ns <- mean(results[,"ns"])
            s <- mean(results[,"s"])
        }
        ## fixed this Dec 7, 2022. I was previously minimizing ONLY on ns.
        if (combiningApproach=="conservativeOld") {
            mostConservativeRow <- which.min(results[,"ns"])[1]
            ns <- results[mostConservativeRow,"ns"]
            s <- results[mostConservativeRow,"s"]
        }
        if (combiningApproach=="conservative") {
            ## when there's a tie based on ns I should ALSO minimize s
            ## take row(s) with minimal ns
            # I previously used which.min - that always takes the FIRST value matching the minimum value
            mostConservativeRows <- which(results[,"ns"] == min(results[,"ns"]))
            results <- results[mostConservativeRows,]
            if(length(mostConservativeRows)>1) {
                ## take row(s) with minimal s
                mostConservativeRows <- which(results[,"s"] == min(results[,"s"]))
                results <- results[mostConservativeRows,]
            }
            # there might still be a tie - I don't care, because all I take is the ns/s counts, which should now have a single value
            results <- unique(results[,c("ns","s")])
            if(dim(results)[1]>1) {
                stop("\n\nERROR - cannot figure out a most conservative solution. code is not working as well as I would expect\n\n")
            }
            ns <- results[,"ns"]
            s <- results[,"s"]
        }
        results <- list(ns=ns, s=s)
        if(extraVerbose) {
            cat("            results AFTER combining:\n")
            print(results)
            cat("\n")
        }
        return(results)
    })
    
    counts_df <- data.frame(codon=1:numCodons,
                            ns=sapply(counts, "[[", "ns"),
                            s=sapply(counts, "[[", "s"))
    return(counts_df)
}    





##### categorizePolymorphisms_new_byCodon - this time we use the sequences directly
# this version now considers one codon at a time
# ancCodons_withUncert - a list, one element per codon, each of which is character vector of possible ancestral codons that we will 'inject' non-ancestral alleles into
# allelesByNuc  - a list, one element per nucleotide, each of which is character vector of alleles
# combiningApproach can be:
#     "mean" (take mean ns and mean s over all possible combinations) or 
#     "conservative" (take the count combination that minimizes ns)

categorizePolymorphisms_new_byCodon <- function(ancCodons_withUncert, 
                                                allelesByNuc, 
                                                combiningApproach="conservative",
                                                extraVerbose=FALSE) {
    if (length(allelesByNuc)/3 != length(ancCodons_withUncert)) {
        stop("ERROR - allelesByNuc and ancCodons_withUncert respresent sequences of different lengths\n\n")    
    }
    if (!combiningApproach %in% c("mean","conservative","conservativeOld") ) {
        stop("\n\nERROR in getCodonChangeCountsFromCodonList - the combiningApproach option can only be mean or conservative or conservativeOld\n\n")
    }
    
    #### work on each codon at a time
    counts <- lapply(1:length(ancCodons_withUncert), function(i) {
        thisCodonStartPos <- 3*(i-1) + 1
        thisCodonEndPos <- thisCodonStartPos+2
        allelesByNucThisCodon <- allelesByNuc[ thisCodonStartPos:thisCodonEndPos ]
        
        ## for each polymorphic position, I will consider it in all possible contexts given polymorphisms at other positions in the same codon, ignoring stop codons
        changesThisCodon <- list()
        # the polymorphic position we are considering is at position j within the codon
        for (j in 1:3) {
            alleles <- allelesByNucThisCodon[[j]]
            if (length(alleles)==1) { # fixed positions
                changesThisCodon[[j]] <- list(ns=0,s=0)
                next
            }
            # remaining positions are polymorphic
            ancCodons <- ancCodons_withUncert[[ i ]]
            otherPositionsSameCodon <- setdiff(1:3, j)
            
            #cat("\ncodon",i,"position",j,"\n")
            #cat("ancCodons",ancCodons,"\n")
            #cat("alleles",alleles,"\n")
            
            ## look at SNPs in the OTHER codon positions to make sure I'm considering all contexts for the current SNP
            # collect all possible pre-polymorphisms codons in allBeforeCodons
            allBeforeCodons <- ancCodons
            for (k in otherPositionsSameCodon) {
                otherNucsSameCodon <- allelesByNucThisCodon[[k]]
                for (ancCodon in ancCodons) {
                    newAncCodon <- ancCodon
                    for (nuc in otherNucsSameCodon) {
                        substring(newAncCodon,k,k) <- nuc
                        if(!newAncCodon %in% allBeforeCodons) { 
                            if(!newAncCodon  %in% stopCodons) {
                                allBeforeCodons <- c(allBeforeCodons, newAncCodon)
                            }
                        }
                    }
                }
                #cat("    otherNucsSameCodon pos ",k,"nucs",otherNucsSameCodon,"\n")
            }
            
            ### then inject each allele of the current polymorphism into each of the allBeforeCodons.
            # we keep any where before and after are different, and don't involve stop codons
            beforeCodonsToCheck <- character()
            afterCodonsToCheck <- character()
            for (allele in alleles) {
                for (beforeCodon in allBeforeCodons) {
                    thisAfterCodon <- beforeCodon
                    substring(thisAfterCodon,j,j) <- allele
                    if (thisAfterCodon %in% stopCodons) { next }
                    if (thisAfterCodon != beforeCodon) { 
                        beforeCodonsToCheck <- c(beforeCodonsToCheck,beforeCodon) 
                        afterCodonsToCheck <- c(afterCodonsToCheck,thisAfterCodon) 
                    }
                }
            }
            if(extraVerbose) {
                cat("beforeCodonsToCheck",beforeCodonsToCheck,"\n")
                cat("afterCodonsToCheck",afterCodonsToCheck,"\n")
            }
            results <- getCodonChangeCounts(beforeCodonsToCheck, afterCodonsToCheck, 
                                            extraVerbose=extraVerbose )
            if(extraVerbose) {
                cat("            results BEFORE combining:\n")
                print(results)
                cat("\n")
            }
            if (combiningApproach=="mean") {
                ns <- mean(results[,"ns"])
                s <- mean(results[,"s"])
            }
            ## fixed this Dec 7, 2022. I was previously minimizing ONLY on ns.
            if (combiningApproach=="conservativeOld") {
                mostConservativeRow <- which.min(results[,"ns"])[1]
                ns <- results[mostConservativeRow,"ns"]
                s <- results[mostConservativeRow,"s"]
            }
            if (combiningApproach=="conservative") {
                ## when there's a tie based on ns I should ALSO minimize s
                ## take row(s) with minimal ns
                mostConservativeRows <- which(results[,"ns"] == min(results[,"ns"]))
                results <- results[mostConservativeRows,]
                if(length(mostConservativeRows)>1) {
                    ## take row(s) with minimal s
                    mostConservativeRows <- which(results[,"s"] == min(results[,"s"]))
                    results <- results[mostConservativeRows,]
                }
                # there might still be a tie - I don't care, because all I take is the ns/s counts, which should now have a single value
                results <- unique(results[,c("ns","s")])
                if(dim(results)[1]>1) {
                    stop("\n\nERROR - cannot figure out a most conservative solution. code is not working as well as I would expect\n\n")
                }
                ns <- results[,"ns"]
                s <- results[,"s"]
            }
            
            results <- list(ns=ns, s=s)
            if(extraVerbose) {
                cat("            results AFTER combining:\n")
                print(results)
                cat("\n")
            }
            changesThisCodon[[j]] <- results
        }
        changesThisCodon_df <- data.frame(Pn=sapply(changesThisCodon, "[[", "ns"), 
                                          Ps=sapply(changesThisCodon, "[[", "s"))
        return(changesThisCodon_df)
    })
    ## add codon index to each of those tables
    counts <- lapply(1:length(counts), function(i) {
        tempDF <- counts[[i]]
        tempDF[,"codon"] <- i
        tempDF[,"codonPos"] <- rownames(tempDF)
        tempDF <- tempDF[,c("codon", "codonPos", "Pn", "Ps")]
        return(tempDF)
    })
    counts_df <- do.call("rbind", counts)
    #print (counts_df)
    #cat("\n\nFinished in categorizePolymorphisms_new_byCodon function\n\n")
    return(counts_df)
}



##### categorizePolymorphisms_new_byNuc - this time we use the sequences directly
# this version still considers one nucleotide at a time, but there are unusual cases where we need to consider all polymorphisms within the same codon at the same time, so I need to change this. An example is in test_unusualPolymorphism.fa
# ancCodons_withUncert - a list, one element per codon, each of which is character vector of possible ancestral codons that we will 'inject' non-ancestral alleles into
# allelesByNuc  - a list, one element per nucleotide, each of which is character vector of alleles
# combiningApproach can be:
#     "mean" (take mean ns and mean s over all possible combinations) or 
#     "conservative" (take the count combination that minimizes ns)
categorizePolymorphisms_new_byNuc <- function(ancCodons_withUncert, 
                                              allelesByNuc, 
                                              combiningApproach="conservative") {
    if (length(allelesByNuc)/3 != length(ancCodons_withUncert)) {
        stop("ERROR - allelesByNuc and ancCodons_withUncert respresent sequences of different lengths\n\n")    
    }
    if (!combiningApproach %in% c("mean","conservative","conservativeOld") ) {
        stop("\n\nERROR in getCodonChangeCountsFromCodonList - the combiningApproach option can only be mean or conservative or conservativeOld\n\n")
    }
    codonsToTake <- rep( 1:length(ancCodons_withUncert), each=3)
    codonPositions <- rep( 1:3, length(ancCodons_withUncert))
    counts <- lapply(1:length(allelesByNuc), function(i) {
        alleles <- allelesByNuc[[i]]
        if (length(alleles)==1) { # fixed positions
            return(list(ns=0,s=0))
        }
        # remaining positions are polymorphic
        ancCodons <- ancCodons_withUncert[[ codonsToTake[i] ]]
        codonPos <- codonPositions[i]
        # arbitrarily take the first allele as the one I will compare the others to, and figure out what the ancestral codon looks like if that allele is present
        firstAllele <- alleles[1]
        ancCodonsFirstAllele <- sapply(ancCodons, function(x, 
                                                           allele1=firstAllele, 
                                                           pos=codonPos) {
            y <- x
            substring(y,pos,pos) <- allele1
            return(y)
        })
        ## look at the other alleles, compare to first
        remainingAlleles <- alleles[2:length(alleles)]
        
        # inject the first allele, and compare each other allele to that, averaging over possible ancestral codons, but summing over alternative alleles (so for triallelic SNPs I should have a total count of 2 polymorphisms)
        countsByAllele <- lapply(remainingAlleles, function(allele) {
            codonsBefore <- character()
            codonsAfter <- character()
            for (ancCodon in ancCodonsFirstAllele) {
                newCodon <- ancCodon
                substring(newCodon, codonPos, codonPos) <- allele
                if(ancCodon == newCodon) {next} # irrelevant to count when no change
                codonsBefore <- c(codonsBefore, ancCodon)
                codonsAfter <- c(codonsAfter,newCodon)
            }
            results <- getCodonChangeCounts(codonsBefore, codonsAfter)
            if (combiningApproach=="mean") {
                ns <- mean(results[,"ns"])
                s <- mean(results[,"s"])
            }

            ## fixed this Dec 7, 2022. I was previously minimizing ONLY on ns.
            if (combiningApproach=="conservativeOld") {
                mostConservativeRow <- which.min(results[,"ns"])[1]
                ns <- results[mostConservativeRow,"ns"]
                s <- results[mostConservativeRow,"s"]
            }
            if (combiningApproach=="conservative") {
                ## when there's a tie based on ns I should ALSO minimize s
                ## take row(s) with minimal ns
                mostConservativeRows <- which(results[,"ns"] == min(results[,"ns"]))
                results <- results[mostConservativeRows,]
                if(length(mostConservativeRows)>1) {
                    ## take row(s) with minimal s
                    mostConservativeRows <- which(results[,"s"] == min(results[,"s"]))
                    results <- results[mostConservativeRows,]
                }
                # there might still be a tie - I don't care, because all I take is the ns/s counts, which should now have a single value
                results <- unique(results[,c("ns","s")])
                if(dim(results)[1]>1) {
                    stop("\n\nERROR - cannot figure out a most conservative solution. code is not working as well as I would expect\n\n")
                }
                ns <- results[,"ns"]
                s <- results[,"s"]
            }

            results <- list(ns=ns, s=s)
            return(results)
        })
        nsThisPosition <- sum(sapply(countsByAllele, "[[", "ns"), na.rm=TRUE)
        sThisPosition <- sum(sapply(countsByAllele, "[[", "s"), na.rm=TRUE)
        output <- list(ns=nsThisPosition, s=sThisPosition)
        return(output)
    })
    counts_df <- data.frame(Pn=sapply(counts, "[[", "ns"), Ps=sapply(counts, "[[", "s"))
    return(counts_df)
}

#### categorizePolymorphisms - figures out whether polymorphisms are synonymous or non
# df is a table like positionTable
categorizePolymorphisms <- function(df, countsTable=aln_tables[["mel"]], species="mel") {
    polyColumn <- paste(species, "Poly", sep="")
    ancColumn <- paste(species, "Anc", sep="")
    if (!polyColumn %in% colnames(df)) {
        stop("ERROR - there should be a column called",polyColumn,"in the dataframe you supplied but there is not\n\n")
    }
    if (!ancColumn %in% colnames(df)) {
        stop("ERROR - there should be a column called",ancColumn,"in the dataframe you supplied but there is not\n\n")
    }
    
    outputCol_N <- paste(species, "_Pn", sep="")
    outputCol_S <- paste(species, "_Ps", sep="")
    df[,outputCol_N] <- NA
    df[,outputCol_S] <- NA
    
    # get each ancestral codon (names are CHARACTERS of the codon position)
    ancNucs <- df[,ancColumn]
    # replace any unresolved ancestral nucs by ? (they have a space in the character)
    ancNucs[grep(" ",ancNucs)] <- "?"
    # now split into codons
    ancCodons <- split(ancNucs,df[,"codon"])
    ancCodons <- lapply(ancCodons, paste, collapse="")
    
    # polyPos: indices of the polymorphic positions
    polyPos <- which(df[,polyColumn])
    
    # alleles = a list, one item for each polymorphic position, contents are the alleles.
    alleles <- lapply(polyPos, function(x) {
        rownames(countsTable)[which(countsTable[,x]>0)]
    })
    
    # codon indices for each polymorphic position
    polyCodonIndices <- df[polyPos,"codon"]
    polyCodonPos <-  df[polyPos,"codon_pos"]  # 1, 2 or 3, for where SNP is in codon
    
    # polyCodonsAncCodons: a list, one item per polymorphism, names are non-unique if there is >1 change in a codon
    polyCodonsAncCodons <- ancCodons[ match(as.character(polyCodonIndices), names(ancCodons)) ]
    
    # for each polymorphic position, take the ancestral codon, inject the FIRST allele, and use that as the allele to compare all the others to
    # (slightly better than using ancestral codon, if the ancestral allele for this nucleotide is ?)
    # and for each other allele, count NS and S, and add them up for this position
    for (i in 1:length(polyPos)) {
        thisPolyPos <- polyPos[i]
        thisPolyCodonPos <- polyCodonPos[i]
        ancAllele <- as.character(df[thisPolyPos,ancColumn])
        ancCodon <- polyCodonsAncCodons[[i]]
        thisPosAlleles <- alleles[[i]]
        
        # get the codon that the first allele encodes
        firstAlleleCodon <- ancCodon
        substring(firstAlleleCodon,thisPolyCodonPos,thisPolyCodonPos) <- thisPosAlleles[1]
        
        #cat("i",i,"polyPos",thisPolyPos,"codonPos",thisPolyCodonPos)
        #cat(" ancAllele",ancAllele,"ancCodon",ancCodon,"\n")
        ns <- 0
        s <- 0
        for(thisAllele in thisPosAlleles[2:length(thisPosAlleles)]) {
            ## this code left over from when I was comparing to ancestral allele not first allele - the following if statement should never be true
            if (thisAllele == thisPosAlleles[1]) {
                # cat("    same as first allele. thisAllele",thisAllele,"\n")
                # next
                stop("    ERROR - check code!\n\n")
            }
            # inject allele into ancestral codon
            newCodon <- firstAlleleCodon
            substring(newCodon,thisPolyCodonPos,thisPolyCodonPos) <- thisAllele
            changes <- getCodonChangeCounts(firstAlleleCodon,newCodon)
            ns <- ns + changes[,"ns"]
            s <- s + changes[,"s"]
            #cat("    derived. thisAllele",thisAllele,"newCodon",newCodon,"\n")
            #cat("      ns",changes[,"ns"],"s",changes[,"s"],"\n")
        }
        df[thisPolyPos,outputCol_N] <- ns
        df[thisPolyPos,outputCol_S] <- s
    }
    return(df)
}

#### for fixed changes: 
# nucDF is a data frame with one row for every position in the alignment
# codonDF is the output of something like getCodonChangeCounts
addFixedCountsToNucPositionTable <- function(nucDF, codonDF, 
                                             outputColPrefix, 
                                             columnsToAdd=c("ns","s"),
                                             outputColBaseNames=c("Dn","Ds") ) {
    if(sum(!columnsToAdd %in% colnames(codonDF))) {
        stop("ERROR - cannot find all these columns in the codon table:",columnsToAdd,"\n\n")
    }
    outputColumns <- paste(outputColPrefix, outputColBaseNames, sep="_")
    for (tempColIndex in 1:length(columnsToAdd)) {
        oldColname <- columnsToAdd[tempColIndex]
        newColname <- outputColumns[tempColIndex]
        rowsToAddTo <-  which(nucDF[,"codon_pos"]==1)
        nucDF[,newColname] <- NA
        nucDF[rowsToAddTo,newColname] <- codonDF[,oldColname]
    }
    return(nucDF)
} 

#### makeCodonsFromSeqAsListWithUncertainty - a utility function to group nucleotides into codons

## mylist is a list, one element per nucleotide position, containing any nucleotide we want to consider at that position
## output is a list, one element per codon, containing any codons that could be made from that position
makeCodonsFromSeqAsListWithUncertainty <- function(mylist,
                                                   excludeStopCodons=stopCodons,
                                                   populationTag=NULL, # to help track errors
                                                   extraVerbose=FALSE) {
    numNT <- length(mylist)
    if (numNT/3 != round(numNT/3)) {
        stop("\n\nERROR - sequence length is not a multiple of three",numNT,"\n")
    }
    numCodons <- numNT/3
    codonStarts <- 1 + 3*((1:numCodons)-1)
    codonEnds <- codonStarts + 3 - 1
    codons <- lapply(1:numCodons, function(i) {
        nucs <- mylist[ codonStarts[i]:codonEnds[i] ]
        allCodonVariants <- character()
        for (nuc1 in nucs[[1]]) {
            for (nuc2 in nucs[[2]]) {
                for (nuc3 in nucs[[3]]) {
                    allCodonVariants <- c(allCodonVariants,paste(nuc1,nuc2,nuc3,sep=""))
                }
            }
        }
        if(!is.null(allCodonVariants)) {
            if(length(setdiff(allCodonVariants, excludeStopCodons))==0) {
                cat("all codons",allCodonVariants,"\n")
                cat("non-stop codons",setdiff(allCodonVariants, excludeStopCodons),"\n")
                errorMsg <- "\n\nERROR"
                if (!is.null(populationTag)) {
                    errorMsg <- paste(errorMsg, " in population ", populationTag, sep="")
                }
                errorMsg <- paste(errorMsg, 
                                  ", in codon ",i,
                                  " - all codon possibilities are stop codons - that's odd\n\n", 
                                  sep="")
                stop(errorMsg)
            }
            allCodonVariants <- setdiff(allCodonVariants, excludeStopCodons)
        }
        # xxx get rid of any that contain N (added this Dec 9 2022)
        # maybe I should be excluding N even before I get to this step
        allCodonVariants <- grep("N", allCodonVariants, invert=TRUE, value=TRUE)
        
        return(allCodonVariants)
    })
    return(codons)
}


##### filterAlnRemoveRareVariants - takes an alignment (from a single population), and for any allele whose frequency is lower than the specified threshold, replaces that allele with the major allele 
# (if there's a tie between major allele freqs we arbitrarily pick one)
filterAlnRemoveRareVariants <- function(myAln, alleleFreqThreshold=0) {
    if(alleleFreqThreshold==0) {
        cat("\n\n    WARNING - you did not specify an allele frequency threshold, so we are not filtering out rare alleles\n\n")
        return(myAln)
    }
    myAln_df <- as.data.frame(as.matrix(myAln), stringsAsFactors=FALSE)
    # frequencies (over all non-N bases)
    myAln_freqs <- getACGTfreqs(tabulateDF(myAln_df))
    # which positions are polymorphic?
    polymorphicPositions <- which(apply(myAln_freqs, 2, function(x) { sum(x>0) > 1 }))
    # go through those, and if any alleles have freq under threshold, modify the alignment
    positionsModified <- 0
    allelesModified <- 0
    for (i in polymorphicPositions) {
        freqs <- myAln_freqs[,i]
        nonZeroFreqs <- freqs[which(freqs>0)]
        # if there are no rare alleles, we don't need to do anything. use <= to mimic MKT website
        if (sum(nonZeroFreqs <= alleleFreqThreshold) == 0) {next}
        # but if there are, we replace those in the alignment with the major allele (and if there is a tie between major alleles we choose arbitrarily)
        positionsModified <- positionsModified + 1
        majorAllele <- rownames(myAln_freqs)[ which.max(freqs)[1]]
        # use <= to mimic MKT website
        rareAlleles <- rownames(myAln_freqs)[ (freqs>0 & freqs<=alleleFreqThreshold)]
        for (rareAllele in rareAlleles) {
            #cat("removing allele",rareAllele,"at position",i,"\n")
            myAln_df[,i] <- gsub(rareAllele,majorAllele,myAln_df[,i])
            allelesModified <- allelesModified + 1
        }
    }
    cat("        modified",allelesModified,"alleles at",positionsModified,"positions\n")
    # turn the data.frame back into a BStringSet
    aln_filt <- BStringSet(apply(myAln_df, 1, function(y) {
        BString(paste(y, collapse=""))
    }))
    return(aln_filt)
}

##### doMKtest is a standalone function, that runs the MK tests given a fasta alignment file and a bunch of options:

### I could:
# (a) run it on a fasta file that I name in the myAlnFile option
# (b) run it on a BSStringSet alignment

# myAlnFile: a fasta file, in-frame alignment
# myAln: a BSStringSet alignment (must also specify outfileStem)

# pop1seqs: character vector, names of sequences in first population (e.g. mel)
# pop2seqs: character vector, names of sequence(s) in second population (e.g. sim)
# polarize: whether to do a polarized test in addition to regular MK test? (if so, must supply outgroups
# outgroupSeqs = list of character vectors, containing outgroup(s) in order I want to use them
# combiningApproach: how to deal with >1 change in the same codon. options are:
#     "mean" (take mean ns and mean s over all possible combinations) or 
#     "conservative" (take the count combination that minimizes ns)
#                     flagRareAlleles=FALSE, alleleFreqThreshold=0
# flagRareAlleles: flags rare polymorphisms with freq <= alleleFreqThreshold, but does not change the MK counts
# filterRareAlleles: remove rare alleles with freq <= alleleFreqThreshold
#        use <= to mimic MKT website
# alleleFreqThreshold: used by flagRareAlleles
# writeAncFasta and writeMKoutput - will I write output files for each individual alignment, or will I wait until I combine output from several alignments (=default) 
doMKtest <- function(myAlnFile=NULL, 
                     myAln=NULL, outfileStem=NULL,
                     outDir=NULL, 
                     pop1seqs=NULL, pop2seqs=NULL, 
                     pop1alias=NULL, pop2alias=NULL,
                     polarize=FALSE, outgroupSeqs=NULL,
                     combiningApproach="conservative", 
                     flagRareAlleles=FALSE, 
                     filterRareAlleles=FALSE, 
                     alleleFreqThreshold=0,
                     writeAncFasta=FALSE, writeMKoutput=FALSE,
                     quiet=FALSE, extraVerbose=FALSE) {
    
    ### some checks
    if(is.null(myAln) & is.null(myAlnFile)) {
        stop("\n\nERROR in doMKtest - must specify either the name of an alignment file, or a BStringSet alignment object\n\n")
    }
    if(!is.null(myAln) & !is.null(myAlnFile)) {
        stop("\n\nERROR in doMKtest - must specify either the name of an alignment file, or a BStringSet alignment object but you have specified BOTH!\n\n")
    }
    if(!is.null(myAlnFile)) {
        if (!file.exists(myAlnFile)) {
            stop("ERROR - file ",myAlnFile," does not exist\n\n")
        }
    }
    if(!is.null(myAln)) {
        if(class(myAln) != "BStringSet") {stop("\n\nERROR - myAln should be a BStringSet object\n\n")}
        if(is.null(outfileStem)) {
            stop("\n\nERROR - when using the doMKtests function on a BStringSet object, must also specify the outfileStem argument\n\n")
        }
    }
    if (!combiningApproach %in% c("mean","conservative","conservativeOld") ) {
        stop("\n\nERROR in getCodonChangeCountsFromCodonList - the combiningApproach option can only be mean or conservative or conservativeOld\n\n")
    }
    
    ### read in alignment (if we're using myAlnFile), figure out output file names
    if(!is.null(myAlnFile)) {
        cat("\n##### reading alignment from file",myAlnFile,"\n")
        outfileStem <- gsub(".fasta$","", myAlnFile)
        outfileStem <- gsub(".fas$","", outfileStem)
        outfileStem <- gsub(".fa$","", outfileStem)
        outfileStem <- paste(outfileStem, sep=".")
        #### read alignment - I keep it as a BStringSet as I might want to use ? character later
        aln <- readBStringSet(myAlnFile)
        # check for a weird situation where we try to read in an existing empty file. It does work but gives a 0-length alignment
        if(length(aln)==0) {
            stop("\n\nERROR - the alignment file you specified exists but is empty: ",myAlnFile,"\n\n")
        }
        # take only the first word of the seq names
        names(aln) <- sapply(strsplit(names(aln)," "),"[[",1)
    }
    if(!is.null(myAln)) { aln <- myAln }
    
    ## check that all seqs in the alignment are the same length
    if (length(unique(width(aln)))>1) {
        stop("\n\nERROR - your alignment looks bad. Sequences are not all the same length as each other\n\n")
    }
    ## make sure everything is upper case
    aln <- BStringSet(toupper(aln))
    
    
    ## check population names don't overlap with each other
    if(length(intersect(pop1seqs, pop2seqs))>0) {
        overlappingNames <- intersect(pop1seqs, pop2seqs)
        overlappingNames <- paste(overlappingNames, collapse=" ")
        stop("\n\nERROR - there are sequence names found in both pop1seqs and pop2seqs - that's not right.\nThe overlapping names are: ",overlappingNames,"\n\n")
    }
    if(!is.null(outgroupSeqs)) {
        outgroupSeqsToTest <- unlist(outgroupSeqs, use.names = FALSE)
        overlappingNames <- intersect(outgroupSeqsToTest, c(pop1seqs, pop2seqs))
        if(length(overlappingNames)>0) {
            overlappingNames <- paste(overlappingNames, collapse=" ")
            stop("\n\nERROR - there are sequence names found in both outgroupSeqs and pop1seqs/pop2seqs - that's not right.\nThe overlapping names are: ",overlappingNames,"\n\n")
        }
    }
    
    if (!is.null(outDir)) {
        if (!dir.exists(outDir)) { 
            # only need to make the output dir if we are going to write either of the output file types
            if(writeAncFasta | writeMKoutput) { dir.create(outDir) }
        }
        # strip off old directory path before adding the new one:
        outfileStem <- strsplit(outfileStem,"/")[[1]]
        outfileStem <- outfileStem[length(outfileStem)]
        outfileStem <- paste(outDir,outfileStem,sep="/")
    }
    outfileAln <- paste(outfileStem, ".plusAncs.fa", sep="") 
    outfileMK <- paste(outfileStem, ".MK.xlsx", sep="")
    
    aln_mat <- as.matrix(aln)
    aln_df <- as.data.frame(aln_mat, stringsAsFactors=FALSE)
    
    ### check I can deal with all characters
    if ( sum(!toupper(unlist(aln_df)) %in% acceptableNucs) > 0) {
        unrecognizedNucs <- unique(unlist(aln_df))
        unrecognizedNucs <- unrecognizedNucs[ which(! unrecognizedNucs %in% acceptableNucs ) ]
        stop("ERROR - there are characters in the alignment I don't know what to do with:", unrecognizedNucs, "\n\n")
    }
    
    #### how many seqs
    numSeqs <- length(aln)
    if (numSeqs < 3) {
        stop("ERROR - not enough seqs for polarized MK test: there are only",numSeqs,"seqs\n\n")
    }
    
    #### alignment length
    numNT <- width(aln)
    ## check for alignments where seqs are different lengths
    if (length(unique(numNT))>1) {
        stop("ERROR - seqs are not all the same length\n\n")
    }
    numNT <- numNT[1]
    ## check for length being multiple of 3
    numCodons <- numNT/3
    if (numCodons != round(numCodons)) {
        stop("ERROR - alignment length is not a multiple of 3: ",numNT,"\n\n")
    }
    
    cat("    splitting alignment into groups and tabulating nucleotides\n")
    
    #### check that all pop1 and pop2 seqs are in the alignment (and possibly outgroups)
    if (sum(!pop1seqs %in% names(aln))>0) {
        missingSeqs <- setdiff(pop1seqs, names(aln))
        cat("ERROR - there are seqs named in the pop1seqs argument that are not in the alignment file.\nFile=",
            myAlnFile,
            "\nmissing seqs:", paste(missingSeqs, collapse=" "),
            "\nseqs present:", paste(names(aln), collapse=" "), "\n")
        stop()
    }
    if (sum(!pop2seqs %in% names(aln))>0) {
        missingSeqs <- setdiff(pop2seqs, names(aln))
        cat("ERROR - there are seqs named in the pop2seqs argument that are not in the alignment file.\nFile=",
            myAlnFile,
            "\nmissing seqs:", paste(missingSeqs, collapse=" "),
            "\nseqs present:", paste(names(aln), collapse=" "), "\n")
        stop()
    }
    if(!is.null(outgroupSeqs)) {
        ## check for zero length
        if(length(unlist(outgroupSeqs))==0) {
            stop("\n\nERROR - you supplied an empty outgroupSeqs object\n\n")
        }
        ## check whether all outgroupSeqs are in the alignment
        if (sum(!unlist(outgroupSeqs) %in% names(aln))>0) {
            missingSeqs <- setdiff(unlist(outgroupSeqs), names(aln))
            cat("ERROR - there are seqs named in the outgroupSeqs argument that are not in the alignment file.\nFile=",
                myAlnFile,
                "\nmissing seqs:", paste(missingSeqs, collapse=" "),
                "\nseqs present:", paste(names(aln), collapse=" "), "\n")
            stop()
        }
    }
    
    #### check no repeated seq names in each argument
    if (length(pop1seqs) != length(unique(pop1seqs))) {
        stop("ERROR - there are one or more repeated names in the pop1seqs argument\n\n")    
    }
    if (length(pop2seqs) != length(unique(pop2seqs))) {
        stop("ERROR - there are one or more repeated names in the pop2seqs argument\n\n")    
    }
    if (length(unlist(outgroupSeqs)) != length(unique(unlist(outgroupSeqs)))) {
        stop("ERROR - there are one or more repeated names in the outgroupSeqs argument\n\n")    
    }
    
    #### reorder alignment
    seqOrder <- c(pop1seqs, pop2seqs)
    if (!is.null(outgroupSeqs)) { seqOrder <- c(seqOrder, unlist(outgroupSeqs)) }
    seqs_not_used <- setdiff( names(aln), seqOrder)
    aln <- aln[seqOrder]
    aln_mat <- as.matrix(aln)
    aln_df <- as.data.frame(aln_mat, stringsAsFactors=FALSE)
    #return(aln) ## xx temp
    #### get each extant codon for whole alignment (so I can get rid of stop codons
    aln_codonsUnique <- alnSlicesUniqueSeqs(aln, sliceWidth=3)
    #return(aln_codonsUnique) ## xx temp
    #return(aln_codonsUnique) ## xx temp!
    codonsHaveOnlyStopsOrNs <- unlist(lapply(aln_codonsUnique, function(x) {
        nonStopOrNcodons <- unique( setdiff( toupper(x), stopCodons ) )
        # also get rid of any codon containing N or -
        nonStopOrNcodons <- grep("N|-", nonStopOrNcodons, invert=TRUE, value=TRUE)
        testOnlyStops <- length(nonStopOrNcodons)==0
        return(testOnlyStops)
    }))
    if(sum(codonsHaveOnlyStopsOrNs)>0) {
        #cat("!!HERE!!\n")
        ## test for internal stop codons:
        internalCodonTests <- codonsHaveOnlyStopsOrNs[1:(length(codonsHaveOnlyStopsOrNs)-1)]
        if(sum(internalCodonTests)>0) {
            failedTests <- which(internalCodonTests)
            stop("\n\nERROR - there are problem codon(s) where there are no non-stop, non-gap, non-N containing codons.\n",
                 "Position(s) of problem codons are: ",failedTests,
                 "\nAll codons found at those positions are ",aln_codonsUnique[failedTests],"\n\n")
        }
        ## the stop codon must be at the end - we can strip it out
        cat("\n\nWARNING - this alignment has a codon at the end with only stop codons - we will strip it out and not count any changes in this codon\n\n")
        numNT <- numNT - 3
        numCodons <- numCodons - 1
        aln <- narrow(aln, start=1, end=numNT)
        aln_mat <- as.matrix(aln)
        aln_df <- as.data.frame(aln_mat, stringsAsFactors=FALSE)
    }
    
    #### split alignment into populations/outgroups
    seqClassesInAlnOrder <- c(rep("pop1", length(pop1seqs)),
                              rep("pop2", length(pop2seqs)))
    if (!is.null(outgroupSeqs)) { 
        seqClassesInAlnOrder <- c(seqClassesInAlnOrder, 
                                  rep("out", length(unlist(outgroupSeqs))) )
        ## one alignment for all outgroups
        outgroups_aln <- aln[unlist(outgroupSeqs)]
        # maybe I have several sets of outgroups - get these objects for each set
        # maybe I don't need all these objects - might be able to tidy up code here
        outgroups_seqClassesInAlnOrder <- rep(1:length(outgroupSeqs), 
                                              sapply(outgroupSeqs,length))
        # list of alignments for each outgroup group
        outgroups_aln_split <- split(outgroups_aln, outgroups_seqClassesInAlnOrder)
    }
    
    ##### allele frequencies (I display these in my by-position output, but do not use them for anything else)
    aln_df_split <- split(aln_df, seqClassesInAlnOrder)
    aln_tables <- lapply(aln_df_split, tabulateDF)
    
    # frequencies (over all non-N bases)
    aln_tables_freq <- lapply( aln_tables, getACGTfreqs )
    
    # for each population, get major and minor alleles and their frequencies, flag alleles below a certain frequency threshold (but I do filtering in a separate step, later)
    cat("    looking at allele frequencies\n")
    aln_tables_majorMinor <- lapply(aln_tables_freq, 
                                    getMinorMajorAlleles, 
                                    flagRareAlleles=flagRareAlleles, 
                                    alleleFreqThreshold=alleleFreqThreshold)
    
    #### split alignment by sequence type
    aln_split <- split(aln, seqClassesInAlnOrder)
    
    #### for each of pop1 and pop2, we replace rare alleles with the major allele
    if(filterRareAlleles) {
        alns_filt <- lapply( c("pop1","pop2"), function(thisPop) {
            cat("    filtering alleles for",thisPop,"\n")
            filterAlnRemoveRareVariants(aln_split[[thisPop]],
                                        alleleFreqThreshold=alleleFreqThreshold)
        })
        names(alns_filt) <- c("pop1","pop2")
        aln_split <- as.list(aln_split)
        ## wierd notation, see my Biostrings github issue https://github.com/Bioconductor/Biostrings/issues/39
        aln_split[["pop1"]] <- alns_filt[["pop1"]]
        aln_split[["pop2"]] <- alns_filt[["pop2"]]
        #return(aln_split)
        aln_split <- BStringSetList(aln_split)
    }
    
    
    #### moving on...
    aln_split["pop1_and_pop2"] <- BStringSetList(c(aln_split[["pop1"]], aln_split[["pop2"]]))
    
    #### get each extant codon in each group
    ## after Dec 9 2022, we remove N-containing codons from the returned codons (I changed alnSlicesUniqueSeqs function)
    aln_split_codonsUnique <- lapply(aln_split, function(x) {
        alnSlicesUniqueSeqs(x)
    })
    if(!is.null(outgroupSeqs)) {
        outgroups_aln_split_CodonsUnique <- lapply(outgroups_aln_split, function(x) {
            alnSlicesUniqueSeqs(x)
        })
    }
    
    #### get each extant nucleotide in each group
    aln_split_PositionsUnique <- lapply(aln_split, function(x) {
        alnSlicesUniqueSeqs(x, sliceWidth=1)
    })
    aln_split_PositionsUnique_withoutNsGaps <- lapply(aln_split_PositionsUnique, function(x) {
        lapply(x, function(y) { y[ which(!y %in% c("N","-"))] } )
    })
    
    if(!is.null(outgroupSeqs)) {
        outgroups_aln_split_PositionsUnique <- lapply(outgroups_aln_split, 
                                                      function(x) {
                                                          alnSlicesUniqueSeqs(x, sliceWidth=1)
                                                      })
    }
    
    #### infer ancestors (including uncertainty)
    # stop codons get excluded from the possibilities in uncertain codons 
    # for now we are treating all outgroups (e.g. yak, sec) as one big group. We polarize pop1 using pop2 if possible, and if it's still unresolved we use the outgroup.
    # but I might want to have a hierarchy of outgroups based on species tree.
    cat("    inferring ancestors\n")
    outgroupsForPop1 <- list(aln_split_PositionsUnique[["pop2"]])
    if (!is.null(outgroupSeqs)) {
        outgroupsForPop1 <- c(outgroupsForPop1, outgroups_aln_split_PositionsUnique)
    }
    outgroupsForPop2 <- list(aln_split_PositionsUnique[["pop1"]])
    if (!is.null(outgroupSeqs)) {
        outgroupsForPop2 <- c(outgroupsForPop2, outgroups_aln_split_PositionsUnique)
    }
    if(extraVerbose) {cat("        pop1\n")}
    pop1_anc_withUncert <- reconstructAncestor_includeUncertainty(
        ingroup=aln_split_PositionsUnique[["pop1"]], 
        outgroups=outgroupsForPop1,
        extraVerbose=extraVerbose )
    pop1_anc_withUncert_codons <- makeCodonsFromSeqAsListWithUncertainty(pop1_anc_withUncert,
                                                                         extraVerbose=extraVerbose, 
                                                                         populationTag="pop1_anc")
    
    if(extraVerbose) {cat("        pop2\n")}
    pop2_anc_withUncert <- reconstructAncestor_includeUncertainty(
        ingroup=aln_split_PositionsUnique[["pop2"]], 
        outgroups=outgroupsForPop2,
        extraVerbose=extraVerbose )
    pop2_anc_withUncert_codons <- makeCodonsFromSeqAsListWithUncertainty(pop2_anc_withUncert,
                                                                         extraVerbose=extraVerbose, 
                                                                         populationTag="pop2_anc")
    
    ### if we are polarizing, get the ancestor to pop1 and pop2 using outgroups
    if (polarize) {
        ##  first make a new ingroup that combines pop1_anc_withUncert and pop2_anc_withUncert 
        pop1_and_pop2_ancs <- lapply(1:numNT, function(i) {
            unique(c(pop1_anc_withUncert[[i]], pop2_anc_withUncert[[i]]))
        })
        if(extraVerbose) {cat("        pop1and2\n")}
        pop1_and_pop2_anc_withUncert <- reconstructAncestor_includeUncertainty(
            ingroup=pop1_and_pop2_ancs, 
            outgroups=list(aln_split_PositionsUnique[["out"]]),
            extraVerbose=extraVerbose)
        pop1_and_pop2_anc_withUncert_codons <- makeCodonsFromSeqAsListWithUncertainty(
            pop1_and_pop2_anc_withUncert,
            extraVerbose=extraVerbose, 
            populationTag="pop1_and_pop2_anc")
        if(extraVerbose) {
            cat("pop1_and_pop2_anc_withUncert:\n")
            print(pop1_and_pop2_anc_withUncert)
            cat("\n")
            cat("pop1_and_pop2_anc_withUncert_codons:\n")
            print(pop1_and_pop2_anc_withUncert_codons)
            cat("\n")
        }
    }
    
    #### add inferred ancestors to the alignment and write it out
    if(writeAncFasta) {
        pop1_anc_forFasta <- getAncForFasta(pop1_anc_withUncert)
        pop2_anc_forFasta <- getAncForFasta(pop2_anc_withUncert)
        ancAln <- BStringSet(c(pop1_anc_forFasta, pop2_anc_forFasta))
        names(ancAln) <- c("pop1_anc","pop2_anc")
        if(polarize) { 
            pop1_and_pop2_anc_forFasta <- getAncForFasta(pop1_and_pop2_anc_withUncert) 
            ancAln <- BStringSet(c(pop1_anc_forFasta, pop2_anc_forFasta, pop1_and_pop2_anc_forFasta))
            names(ancAln) <- c("pop1_anc","pop2_anc","pop1_and_pop2_anc")
        }
        newAln <- c(ancAln,aln)
        writeXStringSet(newAln, outfileAln, format="fasta")
    }
    #### make a table showing what's going on with each alignment position
    cat("    starting output table\n")
    positionTable <- data.frame(pos=1:numNT, 
                                codon=rep(1:numCodons, each=3),
                                codon_pos=rep(1:3,numCodons),
                                pop1_anc=sapply(pop1_anc_withUncert, paste, collapse=" "), 
                                #pop1_poly=sapply(aln_split_PositionsUnique[["pop1"]], length)!=1,
                                pop1_poly=sapply(aln_split_PositionsUnique_withoutNsGaps[["pop1"]], length)!=1,
                                pop1_num_alleles=sapply(aln_split_PositionsUnique[["pop1"]], length),
                                pop2_anc=sapply(pop2_anc_withUncert, paste, collapse=" "),
                                #pop2_poly=sapply(aln_split_PositionsUnique[["pop2"]], length)!=1,
                                pop2_poly=sapply(aln_split_PositionsUnique_withoutNsGaps[["pop2"]], length)!=1,
                                pop2_num_alleles=sapply(aln_split_PositionsUnique[["pop2"]], length))
    positionTable[,"fixed_difference"] <- positionTable[,"pop1_anc"] != positionTable[,"pop2_anc"] # not quite want I want, probably ?
    
    if(polarize) { 
        # might be wrong if the anc includes ? or - ?
        positionTable[,"pop1_and_pop2_anc"] <- sapply(pop1_and_pop2_anc_withUncert, paste, collapse=" ")
    }
    
    # add each outgroup seq individually
    if(!is.null(outgroupSeqs)) {
        for(thisOutgroup in outgroupSeqs) {
            positionTable[,thisOutgroup] <- aln_mat[thisOutgroup,]
        }
        rm(thisOutgroup)
    }
    
    ## adds allele counts and frequencies to positionTable
    positionTable <- addCountsToTable(positionTable, aln_tables[["pop1"]], "pop1")
    positionTable <- addCountsToTable(positionTable, aln_tables_majorMinor[["pop1"]], "pop1", transpose=FALSE)
    
    positionTable <- addCountsToTable(positionTable, aln_tables[["pop2"]], "pop2")
    positionTable <- addCountsToTable(positionTable, aln_tables_majorMinor[["pop2"]], "pop2", transpose=FALSE)
    if(!is.null(outgroupSeqs)) {
        positionTable <- addCountsToTable(positionTable, aln_tables[["out"]], "out")
    }
    
    ## xx here - possibly remove rare alleles from polymorphism count
    ## xx also maybe flag alleles that would be removed and are ancestral - this would make additional apparent fixed changes appear with some ways of counting (although I think I would propose to simply not count the polymorphism, rather than counting it as an additional fixed change)
    
    #### go codon by codon to figure out effect of fixed changes.  polymorphisms I look at individually. fixed changes I consider whole codon at once and use codonPaths. I average over uncertainty in codons
    
    ### unpolarized fixed changes:
    cat("    categorizing population 1 vs 2 fixed changes\n")
    pop1_vs_pop2fixedCounts <- getCodonChangeCountsFromCodonList(
        pop1_anc_withUncert_codons, 
        pop2_anc_withUncert_codons, 
        combiningApproach=combiningApproach, 
        extraVerbose=extraVerbose)

    positionTable <- addFixedCountsToNucPositionTable(positionTable, 
                                                      pop1_vs_pop2fixedCounts,
                                                      outputColPrefix="pop1_vs_pop2")
    
    ### categorizes the (pop1) polymorphisms, adding pop1_Pn and pop1_Ps to positionTable
    cat("    categorizing polymorphisms\n")
    pop1_polyCounts <- categorizePolymorphisms_new_byCodon(pop1_anc_withUncert_codons, 
                                                           aln_split_PositionsUnique[["pop1"]],
                                                           combiningApproach=combiningApproach,
                                                           extraVerbose = extraVerbose)
    colnames(pop1_polyCounts) <- paste("pop1", colnames(pop1_polyCounts), sep="_")
    positionTable <- cbind(positionTable, pop1_polyCounts)
    
    pop2_polyCounts <- categorizePolymorphisms_new_byCodon(pop2_anc_withUncert_codons, 
                                                           aln_split_PositionsUnique[["pop2"]],
                                                           combiningApproach=combiningApproach,
                                                           extraVerbose = extraVerbose)
    colnames(pop2_polyCounts) <- paste("pop2", colnames(pop2_polyCounts), sep="_")
    positionTable <- cbind(positionTable, pop2_polyCounts)
    
    if (polarize) {
        #### polarized, pop1 branch:
        cat("    getting polarized changes\n\n")
        if(extraVerbose) {
            cat("getting pop1_and_pop2_anc_to_pop1_FixedCounts:\n")
        }
        pop1_and_pop2_anc_to_pop1_FixedCounts <- getCodonChangeCountsFromCodonList(
            pop1_and_pop2_anc_withUncert_codons,
            pop1_anc_withUncert_codons, 
            extraVerbose=extraVerbose)
        if(extraVerbose) {
            cat("pop1_and_pop2_anc_to_pop1_FixedCounts:\n")
            print (pop1_and_pop2_anc_to_pop1_FixedCounts)
            cat("\n")
        }
        positionTable <- addFixedCountsToNucPositionTable(positionTable, 
                                                          pop1_and_pop2_anc_to_pop1_FixedCounts,
                                                          "pop1_polarized")
        
        #### polarized, pop2 branch:
        if(extraVerbose) {
            cat("getting pop1_and_pop2_anc_to_pop2_FixedCounts:\n")
        }
        pop1_and_pop2_anc_to_pop2_FixedCounts <- getCodonChangeCountsFromCodonList(
            pop1_and_pop2_anc_withUncert_codons,
            pop2_anc_withUncert_codons, 
            combiningApproach=combiningApproach, 
            extraVerbose=extraVerbose)
        if(extraVerbose) {
            print (pop1_and_pop2_anc_to_pop2_FixedCounts)
            cat("\n")
        }
        positionTable <- addFixedCountsToNucPositionTable(positionTable, 
                                                          pop1_and_pop2_anc_to_pop2_FixedCounts,
                                                          "pop2_polarized")
        
        positionTable[,"total_polarized_Dn"] <- positionTable[,"pop1_polarized_Dn"] + 
            positionTable[,"pop2_polarized_Dn"] 
        positionTable[,"total_polarized_Ds"] <- positionTable[,"pop1_polarized_Ds"] + 
            positionTable[,"pop2_polarized_Ds"] 
    }
    
    cat("    doing MK tests\n")
    MKtable <- MKoutput(positionTable, polarize=polarize)
    
    ## another table that has more info in it, for checking purposes
    if(!quiet) {cat("        seqs_not_used",seqs_not_used,"\n")}
    if (length(seqs_not_used)==0) {seqs_not_used <- ""}
    if (length(seqs_not_used)>1) {seqs_not_used <- paste(seqs_not_used, collapse=",")}
    inputName <- "BStringSet"
    if(!is.null(myAlnFile)) {inputName <- myAlnFile}
    finalOutputTable <- data.frame( input=inputName, 
                                    num_seqs=numSeqs,
                                    num_seqs_pop1=length(pop1seqs),
                                    num_seqs_pop2=length(pop2seqs),
                                    num_seqs_outgroup=length(unlist(outgroupSeqs)),
                                    seqs_not_used=seqs_not_used,
                                    length_NT=numNT,
                                    length_AA=numCodons)
    if(!is.null(outgroupSeqs)) {
        finalOutputTable[,"outgroup"] <- unlist(outgroupSeqs)
    }
    finalOutputTable[,"first_pop1_seq"] <- pop1seqs[1]
    finalOutputTable[,"first_pop2_seq"] <- pop2seqs[1]
    finalOutputTable <- cbind(finalOutputTable, MKtable)
    
    if(writeMKoutput) {
        writeOneMKtestToExcel(finalOutputTable, positionTable, outfileMK,
                              pop1alias=pop1alias, pop2alias=pop2alias) 
    }
    return(list(summary=finalOutputTable, positions=positionTable))
}

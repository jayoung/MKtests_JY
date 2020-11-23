### https://en.wikipedia.org/wiki/McDonald%E2%80%93Kreitman_test
 
require(Biostrings)

## on fast codonPathsFromBioperl.Rdata is in /fh/fast/malik_h/user/jayoung/MKtest/forGithub/
## on mac codonPathsFromBioperl.Rdata is in /Users/jayoung/Desktop/mac_workStuff/mac_MKtests/forGithub
load("scripts/codonPaths/codonPathsFromBioperl.Rdata")

acceptableNucs <- c("A","C","G","T","-","N")

colsForVerticalBorderResultsTables <- c("pop1_vs_pop2_chiSq_pVal", "pop1_vs_pop2_Dn", "pop1_vs_pop2_NI",
                        "polarized_pop1_chiSq_pVal", "polarized_pop1_Dn", "polarized_pop1_NI",
                        "polarized_pop2_chiSq_pVal", "polarized_pop2_Dn", "polarized_pop2_NI")
colsForVerticalBorderPositionTables <- c("pop1_anc", "pop1_A", "pop2_A",
                                         "pop1_vs_pop2_Dn", "pop1_polarized_Dn", "total_polarized_Dn")

### makeCodons is a simple function to take a character vector and return strings of each codon (so output vector will be 1/3 the length)
# thisSeq is a character vector
makeCodons <- function( thisSeq ) {
    codonIndices <- 1:(length(thisSeq)/3)
    codonStarts <- 3*(codonIndices-1)+1
    codonEnds <- codonStarts+2
    codons <- sapply( codonIndices , function(i) {
        paste(thisSeq[ codonStarts[i]:codonEnds[i] ], collapse="")
    })
    return(codons)
}
    
#### alnSlicesUniqueSeqs - a function to take an alignment, split it into codons, and for each codon, return a character vector of the diversity of observed codons at that position
# input = an alignment as BStringSet
# output = a list, one element per codon, character vector of what codons are present
alnSlicesUniqueSeqs <- function(myAln, sliceWidth=3) {
    numSlices <- width(myAln)[1] / sliceWidth
    sliceStarts <- 1 + sliceWidth*((1:numSlices)-1)
    sliceEnds <- sliceStarts + sliceWidth - 1
    slices <- lapply( 1:numSlices, function(i) {
        slice <- subseq(myAln, start=sliceStarts[i], end=sliceEnds[i])
        return(unique(as.character(slice)))
    })
    return(slices)
}


## tabulateDF - if alignment is a table, tabulate the nucleotides at each position
tabulateDF <- function(df) {
    y <- apply(df, 2, function(x) {
        table(factor(x, levels=acceptableNucs))
    })
    y
}

## getACGTfreqs - tabulates frequency of each allele over all ACGT alleles (i.e. ignore N and -). works on the output of tabulateDF
getACGTfreqs <- function(df) {
    apply( df[c("A","C","G","T"),], 2, function(x)  {
        x / sum(x)
    })
}

#### reconstruct ancestors using nucleotide count tables. Useful for polarized tests, but perhaps also just to get mel ancestor using sim as outgroup.
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


####### reconstruct ancestors at individual positions (NOT codons), now including uncertainty
# ingroup is a list like the output of alnSlicesUniqueSeqs for the ingroup we are trying to get ancestor for
# outgroup(s) is a list of lists, for outgroups in order. We use the first outgroup to resolve uncertainty if possible, if not we try the second, etc
reconstructAncestor_includeUncertainty <- function(ingroup, outgroups) {
    #if (length(ingroup[[1]]) != 1) {
    if (nchar(ingroup[[1]][1]) != 1) {
        stop("ERROR - alignment should be split into single nucleotide slices for the reconstructAncestor_includeUncertainty function\n\n")
    }
    numPositions <- length(ingroup)
    ancestor <- lapply(1:numPositions, function(i) {
        ingroupThisPos <- ingroup[[i]]
        if (length(ingroupThisPos)==1) {return(ingroupThisPos)}
        outgroupsThisPos <- lapply(outgroups, "[[", i)
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
    return(ancestor)
}

#### get counts of synon and non-synon changes between sets of codons (character vectors)
getCodonChangeCounts <- function(codonsBefore, codonsAfter, paths=codonPaths) {
    if (length(codonsBefore) != length(codonsAfter)){
        stop("ERROR - codonsBefore and codonsAfter are different lengths\n\n")    
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
    # make sure all codons are in the table
    checkCodons <- unique(codonsJoined[differentCodons])
    if (sum(!checkCodons %in% paths[,"path"])>0) {
        missingCodons <- checkCodons[which(!checkCodons %in% paths[,"path"])]
        stop("ERROR - found a codon missing from the paths table, probably contains a stop codon:",missingCodons,"\n\n")
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

#### getCodonChangeCountsFromCodonList = function to get counts of synon and non-synon changes between sets of codons (character vectors)
# this version allows for ambiguity within a codon, and averages result over that ambiguity
# now, input= two lists to compare. Each list has one element per codon, containing character vector of the codons to compare
# combiningApproach can be:
#     "mean" (take mean ns and mean s over all possible combinations) or 
#     "conservative" (take the count combination that minimizes ns)
getCodonChangeCountsFromCodonList <- function(codonListBefore, codonListAfter, paths=codonPaths, combiningApproach="conservative") {
    if (length(codonListBefore) != length(codonListAfter)){
        stop("\n\nERROR in getCodonChangeCountsFromCodonList - codonListBefore and codonListAfter are different lengths\n\n")    
    }
    if (!combiningApproach %in% c("mean","conservative") ) {
        stop("\n\nERROR in getCodonChangeCountsFromCodonList - the combiningApproach option can only be mean or conservative\n\n")
    }
    numCodons <- length(codonListBefore)
    counts <- lapply(1:numCodons, function(i) {
        before <- codonListBefore[[i]]
        after <- codonListAfter[[i]]
        # generate all possible comparisons between before and after
        numComparisons <- length(before) * length(after)
        beforeAll <- rep( before, length(after))
        afterAll <- rep( after, each=length(before))
        # each possible combination shows up as a row in this result table
        result <- getCodonChangeCounts(beforeAll,afterAll)
        # now get a single ns and s considering all possible combinations
        if (combiningApproach=="mean") {
            ns <- mean(result[,"ns"])
            s <- mean(result[,"s"])
        }
        if (combiningApproach=="conservative") {
            mostConservativeRow <- which.min(result[,"ns"])[1]
            ns <- result[mostConservativeRow,"ns"]
            s <- result[mostConservativeRow,"s"]
        }
        return(list(ns=ns, s=s))
    })
    counts_df <- data.frame(codon=1:numCodons,
                            ns=sapply(counts, "[[", "ns"),
                            s=sapply(counts, "[[", "s"))
    return(counts_df)
}    



##### categorizePolymorphisms_new - this time we use the sequences directly
# ancCodons_withUncert - a list, one element per codon, each of which is character vector of possible ancestral codons that we will 'inject' non-ancestral alleles into
# allelesByNuc  - a list, one element per nucleotide, each of which is character vector of alleles
# combiningApproach can be:
#     "mean" (take mean ns and mean s over all possible combinations) or 
#     "conservative" (take the count combination that minimizes ns)
categorizePolymorphisms_new <- function(ancCodons_withUncert, 
                                        allelesByNuc, 
                                        combiningApproach="conservative") {
    if (length(allelesByNuc)/3 != length(ancCodons_withUncert)) {
        stop("ERROR - allelesByNuc and ancCodons_withUncert respresent sequences of different lengths\n\n")    
    }
    if (!combiningApproach %in% c("mean","conservative") ) {
        stop("\n\nERROR in getCodonChangeCountsFromCodonList - the combiningApproach option can only be mean or conservative\n\n")
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
            if (combiningApproach=="conservative") {
                mostConservativeRow <- which.min(results[,"ns"])[1]
                ns <- results[mostConservativeRow,"ns"]
                s <- results[mostConservativeRow,"s"]
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

##### takes MK results from several alignments and makes a single output file / table
combineMKresults <- function(MKresultList, outFile, outDir=NULL,
                             pop1alias=NULL, pop2alias=NULL,
                             getGeneNames=FALSE, geneNameFile="riniTable2_geneOrder.txt", 
                             rowBordersEachGene=FALSE) {
    
    #### check outDir is present
    if(!is.null(outDir)) {
        if (!dir.exists(outDir)) { dir.create(outDir) }
        outFile <- paste(outDir,outFile,sep="/")
    }
    
    #### combine results into a single table
    results_df <- do.call("rbind", lapply(MKresultList, "[[", "summary"))
    ### some cosmetic stuff - add gene name as its own column
    results_df[,"gene_name"] <- row.names(results_df)
    ## get rini's gene names
    if(getGeneNames) {
        if(!file.exists(geneNameFile)) {
            stop("\n\nERROR - gene name file",geneNameFile,"does not exist\n\n")
        }
        riniGeneOrder <- read.delim(geneNameFile)
        results_df[,"gene_name_Rini"] <- riniGeneOrder[match(results_df[,"gene_name"], riniGeneOrder[,"name3"]), "name1"]
        ## reorder rows
        results_df <- results_df[match(riniGeneOrder[,"name1"], results_df[,"gene_name_Rini"] ),]
        ## reorder columns
        newColOrder <- c("gene_name_Rini", "gene_name")
        newColOrder <- c( newColOrder, setdiff( colnames(results_df) , newColOrder))
        results_df <- results_df[,newColOrder]
    } else {
        ## reorder columns
        newColOrder <- c("gene_name")
        newColOrder <- c( newColOrder, setdiff( colnames(results_df) , newColOrder))
        results_df <- results_df[,newColOrder]
    }

    ####### save output to excel file
    results_df_forOutput <- results_df
    wb <- createWorkbook()
    addDataAndFormatWorksheet(wb, 
                results_df_forOutput, "MKtests", 
                rotateColNames=90, 
                colNamesForVerticalBorder=colsForVerticalBorderResultsTables, 
                pop1alias=pop1alias, pop2alias=pop2alias,
                headerRowHeight=180)
    # rows where new gene occurs have a border above. This is an old thing for when I was looking at Rini's alignments, including sub-alignments for different domains
    if(rowBordersEachGene) {
        tempRowsToTest <- 2:dim(results_df_forOutput)[1]
        tempRowsToAddBorder <- which(results_df_forOutput[tempRowsToTest,"gene name Rini"] != results_df_forOutput[(tempRowsToTest-1),"gene name Rini"])
        tempRowsToAddBorder <- tempRowsToAddBorder + 2
        addStyle(wb, "MKtests", createStyle(border="top"), stack=TRUE, 
                 cols=1:dim(results_df_forOutput)[2], 
                 rows=tempRowsToAddBorder, gridExpand=TRUE) 
    }
    
    #### add positions tables
    positionsTables <- lapply(MKresultList, "[[", "positions")
    for (resultName in names(positionsTables)) {
        posTable <- positionsTables[[resultName]]
        polyCodonsToColor <- posTable[which(posTable[,"pop1_poly"] | 
                                  posTable[,"pop2_poly"] & 
                                  !posTable[,"fixed_difference"]),"codon"]
        fixedCodonsToColor <- posTable[which(posTable[,"fixed_difference"]),"codon"]
        polyPositionsToColor <- which(posTable[,"codon"] %in% polyCodonsToColor)
        fixedPositionsToColor <- which(posTable[,"codon"] %in% fixedCodonsToColor)
        rowIndicesForHorizontalBorder <- which( 
            posTable[2:dim(posTable)[1],"codon"] != 
                posTable[1:(dim(posTable)[1]-1),"codon"] )
        addDataAndFormatWorksheet(wb, 
                    df=posTable, 
                    sheetName=resultName, 
                    rotateColNames=90, numColumnsToFreeze=3,
                    colNamesForVerticalBorder=colsForVerticalBorderPositionTables, 
                    rowIndicesForHorizontalBorder=rowIndicesForHorizontalBorder,
                    rowIndicesColorList=list(green3=polyPositionsToColor, 
                                             red=fixedPositionsToColor),
                    pop1alias=pop1alias, pop2alias=pop2alias,
                    headerRowHeight=150)
    }
    
    #### save Excel file
    saveWorkbook(wb, outFile, overwrite = TRUE)
    
    return(results_df)
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

#### addCountsToTable = utility function to help construct output
addCountsToTable <- function(tableToAddTo, counts, columnPrefix, transpose=TRUE) {
    if(transpose) { x <- t(counts) } else { x <- counts }
    colnames(x) <- paste(columnPrefix, colnames(x), sep="_")
    if( dim(tableToAddTo)[1] != dim(x)[1]) {
        stop("ERROR - tables are not the same size\n\n")
    }
    y <- cbind(tableToAddTo,x)
    return(y)
}

#### makeCodonsFromSeqAsListWithUncertainty - a utility function to group nucleotides into codons

## mylist is a list, one element per nucleotide position, containing any nucleotide we want to consider at that position
## output is a list, one element per codon, containing any codons that could be made from that position
makeCodonsFromSeqAsListWithUncertainty <- function(mylist,
                            excludeStopCodons=c("TAG","TAA","TGA") ) {
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
                stop("\n\nERROR in codon ",i," - all codon possibilities are stop codons - that's odd\n\n")
            }
            allCodonVariants <- setdiff(allCodonVariants, excludeStopCodons)
        }
        return(allCodonVariants)
    })
    codons
}

####### MKoutput works on a table of changes at each nucleodide/codon like positionTable
## output is modelled on Rini's table 2
MKoutput <- function(df, polarize=FALSE, pVal_threshold=0.05) {
    changeTotals <- apply( df[,grep("_[DP][ns]$", colnames(df), value=TRUE, perl=TRUE)], 2, sum, na.rm=TRUE)
    ## numbers for regular MK test (just pop1 versus pop2, no polarization):
    pop1_vs_pop2Nums <- list(Dn=changeTotals["pop1_vs_pop2_Dn"], 
                        Ds=changeTotals["pop1_vs_pop2_Ds"], 
                        Pn=(changeTotals["pop1_Pn"]+changeTotals["pop2_Pn"]), 
                        Ps=(changeTotals["pop1_Ps"]+changeTotals["pop2_Ps"]) )
    output <- data.frame( pop1_vs_pop2_chiSq_pVal=do.call("doChiSq", pop1_vs_pop2Nums),
                pop1_vs_pop2_FET_pVal=do.call("doFisher", pop1_vs_pop2Nums),
                pop1_vs_pop2_Dn=changeTotals["pop1_vs_pop2_Dn"],
                pop1_vs_pop2_Ds=changeTotals["pop1_vs_pop2_Ds"],
                pop1_vs_pop2_Pn=(changeTotals["pop1_Pn"]+changeTotals["pop2_Pn"]),
                pop1_vs_pop2_Ps=(changeTotals["pop1_Ps"]+changeTotals["pop2_Ps"]),
                pop1_vs_pop2_NI=do.call("neutralityIndex",pop1_vs_pop2Nums))
    output[,"pop1_vs_pop2_alpha"] <- 1 - output[,"pop1_vs_pop2_NI"]
    output[,"pop1_vs_pop2_result"] <- testMKresult(output[,"pop1_vs_pop2_FET_pVal"],
                                                   output[,"pop1_vs_pop2_NI"],
                                                   pVal_threshold=pVal_threshold)

    if(polarize) {
        ## numbers for polarized MK test (for now, just changes on pop1 branch versus pop1 polymorphisms, because that is the setup I had for Rini analysis):
        polarized_pop1Nums <- list(Dn=changeTotals["pop1_polarized_Dn"], 
                                  Ds=changeTotals["pop1_polarized_Ds"], 
                                  Pn=changeTotals["pop1_Pn"], 
                                  Ps=changeTotals["pop1_Ps"] )
        polarized_pop2Nums <- list(Dn=changeTotals["pop2_polarized_Dn"], 
                                   Ds=changeTotals["pop2_polarized_Ds"], 
                                   Pn=changeTotals["pop2_Pn"], 
                                   Ps=changeTotals["pop2_Ps"] )
        polarizedOutput <- data.frame( 
                          polarized_pop1_chiSq_pVal=do.call("doChiSq", polarized_pop1Nums),
                          polarized_pop1_FET_pVal=do.call("doFisher", polarized_pop1Nums),
                          polarized_pop1_Dn=changeTotals["pop1_polarized_Dn"], 
                          polarized_pop1_Ds=changeTotals["pop1_polarized_Ds"], 
                          pop1_Pn=changeTotals["pop1_Pn"], 
                          pop1_Ps=changeTotals["pop1_Ps"], 
                          polarized_pop1_NI=do.call("neutralityIndex",polarized_pop1Nums))
        polarizedOutput[,"polarized_pop1_alpha"] <- 1 - polarizedOutput[,"polarized_pop1_NI"]
        polarizedOutput[,"polarized_pop1_result"] <- testMKresult(
            polarizedOutput[,"polarized_pop1_FET_pVal"],
            polarizedOutput[,"polarized_pop1_NI"],
            pVal_threshold=pVal_threshold)

        polarizedOutput[,"polarized_pop2_chiSq_pVal"] <- do.call("doChiSq", polarized_pop2Nums)
        polarizedOutput[,"polarized_pop2_FET_pVal"] <- do.call("doFisher", polarized_pop2Nums)
        polarizedOutput[,"polarized_pop2_Dn"] <- changeTotals["pop2_polarized_Dn"]
        polarizedOutput[,"polarized_pop2_Ds"] <- changeTotals["pop2_polarized_Ds"]
        polarizedOutput[,"pop2_Pn"] <- changeTotals["pop2_Pn"]
        polarizedOutput[,"pop2_Ps"] <- changeTotals["pop2_Ps"]
        polarizedOutput[,"polarized_pop2_NI"] <- do.call("neutralityIndex", 
                                                         polarized_pop2Nums)
        polarizedOutput[,"polarized_pop2_alpha"] <- 1 - polarizedOutput[,"polarized_pop2_NI"]
        polarizedOutput[,"polarized_pop2_result"] <- testMKresult(
            polarizedOutput[,"polarized_pop2_FET_pVal"],
            polarizedOutput[,"polarized_pop2_NI"],
            pVal_threshold=pVal_threshold)
            
        polarizedOutput[,"total_polarized_Dn"] <- changeTotals["total_polarized_Dn"]
        polarizedOutput[,"total_polarized_Ds"] <- changeTotals["total_polarized_Ds"]
        output <- cbind(output, polarizedOutput)
    }
    return(output)
}



##### utility functions for MKoutput
doChiSq <- function(Dn,Ds,Pn,Ps, noWarn=TRUE) {
    if(sum(c(Dn,Ds,Pn,Ps))==0) { return(NA) }
    m <- matrix(c(Dn,Ds,Pn,Ps), nrow=2,ncol=2)
    if (noWarn) {
        return(suppressWarnings(chisq.test(m, correct=FALSE)$p.value))
    } else {
        return(chisq.test(m, correct=FALSE)$p.value)
    }
}
doFisher <- function(Dn,Ds,Pn,Ps) {
    if(sum(c(Dn,Ds,Pn,Ps))==0) { return(NA) }
    m <- matrix(c(Dn,Ds,Pn,Ps), nrow=2,ncol=2)
    fisher.test(m)$p.value
}
neutralityIndex <- function(Dn,Ds,Pn,Ps) {
    NI <- (Pn/Ps) / (Dn/Ds)
    return(NI)
}
alpha <- function(Dn,Ds,Pn,Ps) {
    1 - ((Ds*Pn) / (Dn*Ps))
}
# tests significance
testMKresult <- function(pVal, NI, pVal_threshold) {
    result <- "not signif"
    if (pVal<=pVal_threshold & NI>1) { result <- "negative selection" }
    if (pVal<=pVal_threshold & NI<1) { result <- "positive selection" }
    return(result)
}


addDataAndFormatWorksheet <- function(myWB, df, sheetName, 
                                      rotateColNames=0, headerRowHeight=NULL, 
                                      numColumnsToFreeze=0,
                                      colNamesForVerticalBorder=NULL,
                                      rowIndicesForHorizontalBorder=NULL,
                                      colNamesForBold=NULL, 
    # rowIndicesColorList is a list, where names are colors, and contents are numerical indices of rows to have that color
                                      rowIndicesColorList=NULL, 
    # options above are generic, whereas option below are specific to MK test situation
                                      boldPvalColumns=TRUE,
                                      roundSomeColumns=TRUE,
                                      pop1alias=NULL, pop2alias=NULL) {
    require(openxlsx)
    ## prep the data frame
    df2 <- df
    # in p-value and NI and alpha columns, replace NA/NaN with "N.A."
    for (tempCol in grep("pVal|_NI$|_alpha$",colnames(df2),value=TRUE)) {
        df2[which(is.nan(df2[,tempCol])),tempCol] <- "N.A."
        df2[which(is.na(df2[,tempCol])),tempCol] <- "N.A."
        df2[which(is.infinite(df2[,tempCol])),tempCol] <- "Inf"
    }
    if(boldPvalColumns) {
        pValColNames <- grep("_pVal$", colnames(df2), value=TRUE)
        if(is.null(colNamesForBold)) {
            colNamesForBold <- pValColNames
        } else {
            colNamesForBold <- c(colNamesForBold,pValColNames)
        }
    }
    if(!is.null(colNamesForVerticalBorder)) {
        colIndicesForVerticalBorder <- which(colnames(df2) %in% 
                                                 colNamesForVerticalBorder)
    }
    if (!is.null(colNamesForBold)) {
        colIndicesForBold <- which(colnames(df2) %in% colNamesForBold)
    }

    ## fix up various column names
    colnames(df2) <- gsub("_"," ", colnames(df2))
    ## population aliases
    if (!is.null(pop1alias)) { colnames(df2) <- gsub("pop1", pop1alias, colnames(df2)) } 
    if (!is.null(pop2alias)) { colnames(df2) <- gsub("pop2", pop2alias, colnames(df2)) }
    
    options("openxlsx.borderColour" = "black")
    addWorksheet(myWB, sheetName, zoom=120)
    freezePane(myWB, sheetName, firstRow = TRUE)
    if (numColumnsToFreeze>0) {
        freezePane(myWB, sheetName, 
                   firstActiveRow=2, 
                   firstActiveCol=LETTERS[(numColumnsToFreeze+1)])
    }
    myHeaderStyle <- createStyle(fontName="Arial", fontSize = 12, 
              wrapText=TRUE, textRotation=rotateColNames, 
              textDecoration="bold", valign="top")
    writeData(myWB, sheetName, df2, 
              headerStyle=myHeaderStyle, borderStyle="none")
    ## format header row
    addFilter(myWB, sheetName, row=1, cols=1:ncol(df2))
    if (!is.null(headerRowHeight)) {
        setRowHeights(myWB, sheetName, rows=1, heights=headerRowHeight)
    }
    ## some formatting for ALL cells except header row:
    myStyle <- createStyle(fontName="Arial", fontSize = 12) 
    addStyle(myWB, sheetName, myStyle, stack=TRUE, 
             cols=1:dim(df2)[2], rows=1:dim(df2)[1]+1, gridExpand=TRUE)
    setColWidths(myWB, sheetName, cols=1:dim(df2)[2], widths = "auto")
    # vertical border before selected columns:
    if(!is.null(colNamesForVerticalBorder)) {
        addStyle(myWB, sheetName, createStyle(border="left"), stack=TRUE, 
                 cols=colIndicesForVerticalBorder, 
                 rows=0:dim(df2)[1]+1, gridExpand=TRUE)
    }
    # horizontal border after selected rows:
    if(!is.null(rowIndicesForHorizontalBorder)) {
        addStyle(myWB, sheetName, createStyle(border="bottom"), stack=TRUE, 
                 rows=rowIndicesForHorizontalBorder+1, 
                 cols=1:dim(df2)[2], gridExpand=TRUE)
    }
    # bold for selected columns (each table):
    if(!is.null(colNamesForBold)) {
        addStyle(myWB, sheetName, 
                 createStyle(textDecoration="bold"), stack=TRUE, 
                 cols=colIndicesForBold, 
                 rows=0:dim(df2)[1]+1, gridExpand=TRUE)
    }

    if(!is.null(rowIndicesColorList)) {
        for(thisColor in names(rowIndicesColorList)) {
            addStyle(myWB, sheetName, 
                     createStyle(fontColour=thisColor), stack=TRUE, 
                     rows=(rowIndicesColorList[[thisColor]]+1), 
                     cols=0:dim(df2)[2]+1, gridExpand=TRUE)
        }
    }

    # round some columns. doesn't work right if there are some text values (e.g NA or Inf). There is some sort of answer to this on the internet (search for something like 'openxls format numbers stored as text') but it is too convoluted to be worth the effort
    if(roundSomeColumns) {
        columnsToRound <- grep(" pVal| NI| alpha", colnames(df2))
        if (length(columnsToRound)>0) {
            addStyle(myWB, sheetName, createStyle(numFmt="0.000"), stack=TRUE, 
                     cols=columnsToRound, 
                     rows=0:dim(df2)[1]+1, 
                     gridExpand=TRUE)
        }
    }
}

## getAncForFasta = a small function to take a LIST (one position each element, contains possible uncertain positions) and replace uncertain positions with "?" and make a string of the DNA sequence
getAncForFasta <- function(myanc) {
    dna <- myanc
    dna[which(sapply(myanc, length)>1)] <- "?"
    dna <- unlist(dna)
    dna <- paste(dna,collapse="")
    return(dna)
}

## getMinorMajorAlleles - from a table of allele freqs (only ACGT) report major allele(s) and minor allele(s) and their frequencies. If there is an exact tie there can be >1 major allele

getMinorMajorAlleles <- function(df, flagRareAlleles=FALSE, alleleFreqThreshold=0) {
    bases <- rownames(df)
    output <- apply(df, 2, function(x) {
        majorAlleles <- bases[ which(x==max(x)) ]
        major_allele_freqs <- x[majorAlleles]
        presentAlleles <- bases[ which(x>0) ]
        minorAlleles <- setdiff(presentAlleles, majorAlleles)
        if(length(minorAlleles)==0) { 
            minorAlleles <- NA
            minor_allele_freqs <- NA
            if (flagRareAlleles) {minorAlleleFlags <- NA}
        } else {
            minor_allele_freqs <- x[minorAlleles]
            if (flagRareAlleles) {
                minorAlleleFlags <- rep("common", length(minor_allele_freqs))
                minorAlleleFlags[which(minor_allele_freqs <= alleleFreqThreshold)] <- "rare"
            }
            # round AFTER testing frequency threshold
            minor_allele_freqs <- round(minor_allele_freqs,3)
            if(length(minorAlleles)>1) {
                minorAlleles <- paste(minorAlleles, collapse=" ")
                if (flagRareAlleles) {minorAlleleFlags <- paste(minorAlleleFlags, collapse=" ")}
                minor_allele_freqs <- paste(minor_allele_freqs, collapse=" ")
            }
        }
        if(length(majorAlleles)>0) {
            majorAlleles <- paste(majorAlleles, collapse=" ")
            major_allele_freqs <- paste(round(major_allele_freqs,3), collapse=" ")
        } else {  ## sometimes all seqs in this group have gap
            majorAlleles <- "-"
            major_allele_freqs <- 1
        }
        eachColumnOutput <- list(major=majorAlleles, minor=minorAlleles, 
                majorFreqs=major_allele_freqs, minorFreqs=minor_allele_freqs)
        if (flagRareAlleles) {
            eachColumnOutput[["minorAlleleFlags"]] <- minorAlleleFlags
        }
        return(eachColumnOutput)
    })
    # return(output) ## xx temp
    ## output for now may be ugly. perhaps deal with it later
    output_table <- data.frame(major_allele=sapply(output, "[[", "major"),
                         minor_alleles=sapply(output, "[[", "minor"), 
                         major_allele_freqs=sapply(output, "[[", "majorFreqs"),
                         minor_allele_freqs=sapply(output, "[[", "minorFreqs"))
    if (flagRareAlleles) {
        output_table[,"rareFlags"] <- sapply(output, "[[", "minorAlleleFlags")
    }
    return(output_table)
}




##### doMKtest is a standalone function, that runs the MK tests given a fasta alignment file and a bunch of options:
# myAlnFile: a fasta file, in-frame alignment
# pop1seqs: character vector, names of sequences in first population (e.g. mel)
# pop2seqs: character vector, names of sequence(s) in second population (e.g. sim)
# polarize: whether to do a polarized test in addition to regular MK test? (if so, must supply outgroups
# outgroupSeqs = list of character vectors, containing outgroup(s) in order I want to use them
# combiningApproach: how to deal with >1 change in the same codon. options are:
#     "mean" (take mean ns and mean s over all possible combinations) or 
#     "conservative" (take the count combination that minimizes ns)
#                     flagRareAlleles=FALSE, alleleFreqThreshold=0
# flagRareAlleles: for now this only FLAGS rare polymorphisms - later I want to actually deal with them
# alleleFreqThreshold: used by flagRareAlleles

doMKtest <- function(myAlnFile, outDir=NULL,
                     pop1seqs=NULL, pop2seqs=NULL, 
                     pop1alias=NULL, pop2alias=NULL,
                     polarize=FALSE, outgroupSeqs=NULL,
                     combiningApproach="conservative", 
                     flagRareAlleles=FALSE, alleleFreqThreshold=0,
                     writeAncFasta=FALSE) {
    
    require(openxlsx)
    ### some checks
    if (!file.exists(myAlnFile)) {
        stop("ERROR - file ",myAlnFile," does not exist\n\n")
    }
    if (!combiningApproach %in% c("mean","conservative") ) {
        stop("\n\nERROR in getCodonChangeCountsFromCodonList - the combiningApproach option can only be mean or conservative\n\n")
    }
    
    ### read in alignment, figure out output file names
    cat("\n##### reading alignment from file",myAlnFile,"\n")
    outfileStem <- gsub(".fasta$","", myAlnFile)
    outfileStem <- gsub(".fas$","", outfileStem)
    outfileStem <- gsub(".fa$","", outfileStem)
    outfileStem <- paste(outfileStem, sep=".")
    if (!is.null(outDir)) {
        if (!dir.exists(outDir)) { dir.create(outDir) }
        # strip off old directory path before adding the new one:
        outfileStem <- strsplit(outfileStem,"/")[[1]]
        outfileStem <- outfileStem[length(outfileStem)]
        outfileStem <- paste(outDir,outfileStem,sep="/")
    }
    outfileAln <- paste(outfileStem, ".plusAncs.fa", sep="") 
    outfileMK <- paste(outfileStem, ".MK.xlsx", sep="")
    
    #### read alignment - I keep it as a BStringSet as I might want to use ? character later
    aln <- readBStringSet(myAlnFile)
    # take only the first word of the seq names
    names(aln) <- sapply(strsplit(names(aln)," "),"[[",1)
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
    
    #### get each extant codon for whole alignment (so I can get rid of stop codons
    aln_codonsUnique <- alnSlicesUniqueSeqs(aln, sliceWidth=3)
    
    codonsHaveOnlyStopsOrNs <- unlist(lapply(aln_codonsUnique, function(x) {
        nonStopOrNcodons <- unique( setdiff( toupper(x), c("TAA","TAG","TGA") ) )
        # also get rid of any codon containing N
        nonStopOrNcodons <- grep("N", nonStopOrNcodons, invert=TRUE, value=TRUE)
        testOnlyStops <- length(nonStopOrNcodons)==0
        return(testOnlyStops)
    }))
    
    if(sum(codonsHaveOnlyStopsOrNs)>0) {
        ## test for internal stop codons:
        internalCodonTests <- codonsHaveOnlyStopsOrNs[1:(length(codonsHaveOnlyStopsOrNs)-1)]
        if(sum(internalCodonTests)>0) {
            stop("\n\nERROR - there are internal stop codons in file",myAlnFile,"\n")
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
        # maybe I don't need all these things
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
    
    # get major and minor alleles and their frequencies, flag alleles below a certain frequency threshold
    cat("    looking at allele frequencies\n")
    aln_tables_majorMinor <- lapply(aln_tables_freq, getMinorMajorAlleles, flagRareAlleles=flagRareAlleles, alleleFreqThreshold=alleleFreqThreshold)
    cat("        done looking at allele frequencies\n")
   
    #### split alignment by sequence type
    aln_split <- split(aln, seqClassesInAlnOrder)
    aln_split["pop1_and_pop2"] <- BStringSetList(c(aln_split[["pop1"]], aln_split[["pop2"]]))
    
    #### get each extant codon in each group
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
    #cat("outgroupsForPop1",paste(outgroupsForPop1,sep=" "),"\n")
    #cat("outgroupsForPop2",paste(outgroupsForPop2,sep=" "),"\n")
    #### something is wrong here when I polarize?
    
    #return(list(ingroup=aln_split_PositionsUnique[["pop1"]], 
    #            outgroups=outgroupsForPop1))
    pop1_anc_withUncert <- reconstructAncestor_includeUncertainty(
        ingroup=aln_split_PositionsUnique[["pop1"]], 
        outgroups=outgroupsForPop1 )
    pop1_anc_withUncert_codons <- makeCodonsFromSeqAsListWithUncertainty(pop1_anc_withUncert)
    
    pop2_anc_withUncert <- reconstructAncestor_includeUncertainty(
        ingroup=aln_split_PositionsUnique[["pop2"]], 
        outgroups=outgroupsForPop2 )
    pop2_anc_withUncert_codons <- makeCodonsFromSeqAsListWithUncertainty(pop2_anc_withUncert)
    
    ### if we are polarizing, get the ancestor to pop1 and pop2 using outgroups
    if (polarize) {
        ##  first make a new ingroup that combines pop1_anc_withUncert and pop2_anc_withUncert 
        pop1_and_pop2_ancs <- lapply(1:numNT, function(i) {
            unique(c(pop1_anc_withUncert[[i]], pop2_anc_withUncert[[i]]))
        })
        pop1_and_pop2_anc_withUncert <- reconstructAncestor_includeUncertainty(
            ingroup=pop1_and_pop2_ancs, 
            outgroups=list(aln_split_PositionsUnique[["out"]]) )
        pop1_and_pop2_anc_withUncert_codons <- makeCodonsFromSeqAsListWithUncertainty(pop1_and_pop2_anc_withUncert)
    }
    
    #### add inferred ancestors to the alignment and write it out
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
    
    #### make a table showing what's going on with each alignment position
    cat("    starting output table\n")
    positionTable <- data.frame(pos=1:numNT, 
                codon=rep(1:numCodons, each=3),
                codon_pos=rep(1:3,numCodons),
                pop1_anc=sapply(pop1_anc_withUncert, paste, collapse=" "), 
                pop1_poly=sapply(aln_split_PositionsUnique[["pop1"]], length)!=1,
                pop1_num_alleles=sapply(aln_split_PositionsUnique[["pop1"]], length),
                pop2_anc=sapply(pop2_anc_withUncert, paste, collapse=" "),
                pop2_poly=sapply(aln_split_PositionsUnique[["pop2"]], length)!=1,
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
    cat("    categorizing population 1 vs 2 fixed changes\n")
    ### unpolarized fixed changes:
    pop1_vs_pop2fixedCounts <- getCodonChangeCountsFromCodonList(
        pop1_anc_withUncert_codons, 
        pop2_anc_withUncert_codons, 
        combiningApproach=combiningApproach)
    positionTable <- addFixedCountsToNucPositionTable(positionTable, 
                                                      pop1_vs_pop2fixedCounts,
                                                      outputColPrefix="pop1_vs_pop2")
    
    ### categorizes the (pop1) polymorphisms, adding pop1_Pn and pop1_Ps to positionTable
    cat("    categorizing polymorphisms\n")
    pop1_polyCounts <- categorizePolymorphisms_new(pop1_anc_withUncert_codons, 
                                                  aln_split_PositionsUnique[["pop1"]],
                                                  combiningApproach=combiningApproach)
    colnames(pop1_polyCounts) <- paste("pop1", colnames(pop1_polyCounts), sep="_")
    positionTable <- cbind(positionTable, pop1_polyCounts)
    
    pop2_polyCounts <- categorizePolymorphisms_new(pop2_anc_withUncert_codons, 
                                                  aln_split_PositionsUnique[["pop2"]],
                                                  combiningApproach=combiningApproach)
    colnames(pop2_polyCounts) <- paste("pop2", colnames(pop2_polyCounts), sep="_")
    positionTable <- cbind(positionTable, pop2_polyCounts)
    
    if (polarize) {
        #### polarized, pop1 branch:
        cat("    getting polarized changes\n")
        pop1_and_pop2_anc_to_pop1_FixedCounts <- getCodonChangeCountsFromCodonList(
            pop1_and_pop2_anc_withUncert_codons,
            pop1_anc_withUncert_codons, 
            combiningApproach=combiningApproach)
        positionTable <- addFixedCountsToNucPositionTable(positionTable, 
                                                          pop1_and_pop2_anc_to_pop1_FixedCounts,
                                                          "pop1_polarized")
        
        #### polarized, pop2 branch:
        pop1_and_pop2_anc_to_pop2_FixedCounts <- getCodonChangeCountsFromCodonList(
            pop1_and_pop2_anc_withUncert_codons,
            pop2_anc_withUncert_codons, 
            combiningApproach=combiningApproach)
        positionTable <- addFixedCountsToNucPositionTable(positionTable, 
                                                          pop1_and_pop2_anc_to_pop2_FixedCounts,
                                                          "pop2_polarized")
        
        positionTable[,"total_polarized_Dn"] <- positionTable[,"pop1_polarized_Dn"] + 
            positionTable[,"pop2_polarized_Dn"] 
        positionTable[,"total_polarized_Ds"] <- positionTable[,"pop1_polarized_Ds"] + 
            positionTable[,"pop2_polarized_Ds"] 
    }
    
    #return(positionTable)
    cat("    doing MK tests\n")
    MKtable <- MKoutput(positionTable, polarize=polarize)
    
    ## another table that has more info in it, for checking purposes
    cat("        seqs_not_used",seqs_not_used,"\n")
    if (length(seqs_not_used)==0) {seqs_not_used <- ""}
    if (length(seqs_not_used)>1) {seqs_not_used <- paste(seqs_not_used, collapse=",")}
    finalOutputTable <- data.frame( input=myAlnFile, 
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
    
    ### use openxlsx to write an excel file for MKtable and positionTable (two sheets)
    cat("    saving tables to Excel file\n")
    finalOutputTableForExcel <- finalOutputTable
    positionTableForExcel <- positionTable
    polyCodonsToColor <- positionTable[which(positionTable[,"pop1_poly"] | 
                                            positionTable[,"pop2_poly"] & 
                                            !positionTable[,"fixed_difference"]),"codon"]
    fixedCodonsToColor <- positionTable[which(positionTable[,"fixed_difference"]),"codon"]
    polyPositionsToColor <- which(positionTable[,"codon"] %in% polyCodonsToColor)
    fixedPositionsToColor <- which(positionTable[,"codon"] %in% fixedCodonsToColor)
    rowIndicesForHorizontalBorder <- which( 
        positionTable[2:dim(positionTable)[1],"codon"] !=
        positionTable[1:(dim(positionTable)[1]-1),"codon"] )
    ### make the excel sheet
    wb <- createWorkbook()
    addDataAndFormatWorksheet(wb, 
                finalOutputTableForExcel, "MKtable", 
                rotateColNames=90, 
                colNamesForVerticalBorder=colsForVerticalBorderResultsTables,
                pop1alias=pop1alias, pop2alias=pop2alias,
                headerRowHeight=180)
    addDataAndFormatWorksheet(wb, 
                positionTableForExcel, "nucleotidePositions", 
                rotateColNames=90, numColumnsToFreeze=3,
                colNamesForVerticalBorder=colsForVerticalBorderPositionTables, 
                rowIndicesForHorizontalBorder=rowIndicesForHorizontalBorder, 
                rowIndicesColorList=list(green3=polyPositionsToColor, 
                                         red=fixedPositionsToColor),
                pop1alias=pop1alias, pop2alias=pop2alias,
                headerRowHeight=150)
    saveWorkbook(wb, outfileMK, overwrite = TRUE) # , keepNA = TRUE
    rm(wb)
    return(list(summary=finalOutputTable, positions=positionTable))
}


prepMKdataForPositionPlot <- function(df) {
    ## sum polymorphism counts by codon
    pop1_Pn_by_codon <- tapply(df[,"pop1_Pn"], df[,"codon"], sum)
    pop1_Ps_by_codon <- tapply(df[,"pop1_Ps"], df[,"codon"], sum)
    pop2_Pn_by_codon <- tapply(df[,"pop2_Pn"], df[,"codon"], sum)
    pop2_Ps_by_codon <- tapply(df[,"pop2_Ps"], df[,"codon"], sum) 
    df_sumPolyByCodon <- data.frame(
        codon=as.integer(names(pop1_Pn_by_codon)),
        pop1_Pn=as.integer(pop1_Pn_by_codon),
        pop1_Ps=as.integer(pop1_Ps_by_codon),
        pop2_Pn=as.integer(pop2_Pn_by_codon),
        pop2_Ps=as.integer(pop2_Ps_by_codon)
    )
    df_sumPolyByCodon[,"bothPops_Pn"] <- df_sumPolyByCodon[,"pop1_Pn"] + 
        df_sumPolyByCodon[,"pop2_Pn"]
    df_sumPolyByCodon[,"bothPops_Ps"] <- df_sumPolyByCodon[,"pop1_Ps"] + 
        df_sumPolyByCodon[,"pop2_Ps"]
    
    # fixed change counts are already by codon
    df_fixedByCodon <- df[which(df[,"codon_pos"]==1),]
    df_fixedByCodon <- df_fixedByCodon[, grep("_D", colnames(df_fixedByCodon))]
    
    new_df <- cbind(df_sumPolyByCodon,df_fixedByCodon)
    return(new_df)
}

makeOneMKplot <- function(df, 
                          colNames=list(Dn="pop1_vs_pop2_Dn",
                                        Ds="pop1_vs_pop2_Ds",
                                        Pn="bothPops_Pn",
                                        Ps="bothPops_Ps"),
                          myTitle="unpolarized",
                          myColors=list(Dn="red2", Ds="blue2", 
                                        Pn="grey50", Ps="grey90"), 
                          myLwd=2) {
    ## determine x and y limits
    seqLen <- dim(df)[1]
    yMax <- max( c(df[, colNames[["Dn"]] ], df[, colNames[["Ds"]] ],
                   df[, colNames[["Pn"]] ], df[, colNames[["Ps"]] ]), na.rm=TRUE)
    yMax <- yMax + 1
    
    ## blank plot, with a horizontal line for the protein
    plot(x=c(1,seqLen*1.1), y=c((0-yMax),yMax), "n", 
         main=myTitle, xlab="pos (AA)", ylab="num changes per codon", 
         bty="n", xaxt="n", las=2, mgp=c(3,0.75,0))
    segments(x0=1, y0=0, x1=seqLen, y1=0)
    ## prot length at right
    text(x=(seqLen*1.02), y=0, labels=paste(seqLen,"aa"), adj=0, cex=0.75)
    ## label at y=halfmax for fixed/poly
    mtext("fixed", side=2, at=yMax/2, adj=0.5, line=1.75, srt=90, cex=0.75)
    mtext("poly", side=2, at=(0-yMax/2), adj=0.5, line=1.75, srt=90, cex=0.75)
    ## non-synon fixed
    if (sum(df[, colNames[["Dn"]] ]>0)) {
        Dn_table <- df[,c("codon",colNames[["Dn"]])] 
        Dn_table <- Dn_table[which(Dn_table[,colNames[["Dn"]] ]>0),]
        segments(x0=Dn_table[,"codon"], y0=0,
                 x1=Dn_table[,"codon"], y1=Dn_table[,colNames[["Dn"]] ],
                 col=myColors[["Dn"]], lwd=myLwd)
    }
    ## synon fixed 
    if (sum(df[,colNames[["Ds"]] ]>0)) {
        Ds_table <- df[,c("codon",colNames[["Dn"]],colNames[["Ds"]])] 
        Ds_table <- Ds_table[which(Ds_table[,colNames[["Ds"]] ]>0),]
        # for y1 we start at pop1_vs_pop2_Dn, so that values will be stacked for codons that have both Dn and Ds
        segments(x0=Ds_table[,"codon"], y0=Ds_table[,colNames[["Dn"]]],
                 x1=Ds_table[,"codon"], y1=Ds_table[,colNames[["Ds"]] ],
                 col=myColors[["Ds"]], lwd=myLwd)
    }
    ## non-synon poly
    if (sum(df[,colNames[["Pn"]]]>0)) {
        Pn_table <- df[which(df[, colNames[["Pn"]] ]>0),]
        segments(x0=Pn_table[,"codon"], y0=0,
                 x1=Pn_table[,"codon"], y1=(0-Pn_table[, colNames[["Pn"]] ]),
                 col=myColors[["Pn"]], lwd=myLwd)
    }
    ## synon poly
    if (sum(df[, colNames[["Ps"]] ]>0)) {
        Ps_table <- df[which(df[, colNames[["Ps"]] ]>0),]
        segments(x0=Ps_table[,"codon"], y0=(0-Ps_table[,colNames[["Pn"]] ]),
                 x1=Ps_table[,"codon"], y1=(0-Ps_table[,colNames[["Ps"]] ]),
                 col=myColors[["Ps"]], lwd=myLwd)
    }
    
    ## legends
    legend("topright", legend=c("Dn","Ds"), 
           col=c(myColors[["Dn"]], myColors[["Ds"]]), 
           lwd=myLwd, cex=0.75)
    legend("bottomright", legend=c("Pn","Ps"), 
           col=c(myColors[["Pn"]], myColors[["Ps"]]),
           lwd=myLwd, cex=0.75)    
}

## df is a position table output by doMKtest
plotMKpositions <- function(df, title=NULL,
                            plotPolarizedPop1=FALSE, 
                            plotPolarizedPop2=FALSE, 
                            setNumPlots=TRUE, ## will want this FALSE if we set mfrow outside the function
                            pop1alias=NULL, pop2alias=NULL,
                            myColors=list(Dn="red2", Ds="blue2", 
                                          Pn="grey50", Ps="grey90")) {
    
    if(setNumPlots) {
        numPlots <- 1
        if(plotPolarizedPop1) { numPlots <- numPlots + 1 }
        if(plotPolarizedPop2) { numPlots <- numPlots + 1 }
        par(mfrow=c(numPlots, 1))
    }
    par(mar=c(2.1, 4.1, 2.1, 2.1))
    
    ## get data in a better format
    new_df <- prepMKdataForPositionPlot(df)
    
    ## do the unpolarized plot
    titleUnpolarized <- "unpolarized"
    if(!is.null(title)) { titleUnpolarized <- paste(title,titleUnpolarized, sep=" - ") }
    makeOneMKplot(new_df, myTitle=titleUnpolarized)
    
    ## polarized plot for population 1
    if(plotPolarizedPop1) { 
        titlePop1 <- "polarized pop1"
        if (!is.null(title)) { titlePop1 <- paste(title,titlePop1, sep=" - ") }
        if (!is.null(pop1alias)) { titlePop1 <- gsub("pop1", pop1alias, titlePop1) } 
        makeOneMKplot(new_df, 
                      colNames=list(Dn="pop1_polarized_Dn",
                                    Ds="pop1_polarized_Ds",
                                    Pn="pop1_Pn",
                                    Ps="pop1_Ps"),
                      myTitle=titlePop1)
    }
    
    ## polarized plot for population 2
    if(plotPolarizedPop2) { 
        titlePop2 <- "polarized pop2"
        if(!is.null(title)) { titlePop2 <- paste(title,titlePop2, sep=" - ") }
        if (!is.null(pop2alias)) { titlePop2 <- gsub("pop2", pop2alias, titlePop2) } 
        makeOneMKplot(new_df, 
                      colNames=list(Dn="pop2_polarized_Dn",
                                    Ds="pop2_polarized_Ds",
                                    Pn="pop2_Pn",
                                    Ps="pop2_Ps"),
                      myTitle=titlePop2)
    }
    return(NULL)
}


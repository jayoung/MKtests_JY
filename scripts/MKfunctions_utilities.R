
acceptableNucs <- c("A","C","G","T","-","N")
stopCodons <- c("TAA","TAG","TGA")

colsForVerticalBorderResultsTables <- c("pop1_vs_pop2_chiSq_pVal", "pop1_vs_pop2_Dn", "pop1_vs_pop2_NI",
                                        "polarized_pop1_chiSq_pVal", "polarized_pop1_Dn", "polarized_pop1_NI",
                                        "polarized_pop2_chiSq_pVal", "polarized_pop2_Dn", "polarized_pop2_NI")
colsForVerticalBorderPositionTables <- c("pop1_anc", "pop1_A", "pop2_A",
                                         "pop1_vs_pop2_Dn", "pop1_polarized_Dn", "total_polarized_Dn")


### makeCodons is a simple function to take a character vector, where each element is a single nucleotide, and return strings of each codon (so output vector will be 1/3 the length)
# thisSeq is a character vector
makeCodons <- function( thisSeq ) {
    
    ### some checks:
    if(class(thisSeq) != "character") {stop("\n\nERROR - makeCodons takes only character vectors as input\n\n")}
    if (unique(nchar(thisSeq))[1] != 1) {
        stop("\n\nthe items in the input sequences should each have one character in them - you may need to split the seq into individual nucleotides before using makeCodons\n\n")
    }
    if (length(unique(nchar(thisSeq)))>1) {
        stop("\n\nthe items in the input sequences should each have one character in them - you may need to split the seq into individual nucleotides before using makeCodons\n\n")
    }
    if(length(thisSeq)<3) {stop("\n\nERROR - there's <3 nucleotides here, cannot make any codons\n\n")}
    
    ### actually do it
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


## tabulateDF - for alignments that I've converted to a table, tabulate the nucleotides at each position
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
    if(is.na(pVal)) {return(NA)}
    if(is.na(pVal)) {return(NA)}
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
getMinorMajorAlleles <- function(df, 
                                 flagRareAlleles=FALSE, 
                                 alleleFreqThreshold=0) {
    bases <- rownames(df)
    # look at each column of the data.frame, i.e. each alignment position
    output <- apply(df, 2, function(x) {
        #cat("position",x,"\n")
        majorAlleles <- bases[ which(x==max(x)) ]
        major_allele_freqs <- x[majorAlleles]
        presentAlleles <- bases[ which(x>0) ]
        minorAlleles <- setdiff(presentAlleles, majorAlleles)
        # initialize rareAlleles - will populate it later
        rareAlleles <- c()
        if(length(minorAlleles)==0) { # there was no polymorphism
            minorAlleles <- NA
            minor_allele_freqs <- NA
            if (flagRareAlleles) {minorAlleleFlags <- NA}
        } else { # there is polymorphism
            minor_allele_freqs <- x[minorAlleles]
            if (flagRareAlleles) {
                minorAlleleFlags <- rep("common", length(minor_allele_freqs))
                # use <= to mimic MKT website
                rareAlleleTest <- which(minor_allele_freqs <= alleleFreqThreshold)
                if(length(rareAlleleTest)>0) {
                    minorAlleleFlags[rareAlleleTest] <- "rare"
                }
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
        eachColumnOutput <- list(major=majorAlleles, 
                                 minor=minorAlleles, 
                                 majorFreqs=major_allele_freqs, 
                                 minorFreqs=minor_allele_freqs)
        if (flagRareAlleles) {
            eachColumnOutput[["minorAlleleFlags"]] <- minorAlleleFlags
        }
        return(eachColumnOutput)
    })
    ## turn it into a table 
    output_table <- data.frame(major_allele=sapply(output, "[[", "major"),
                               minor_alleles=sapply(output, "[[", "minor"), 
                               major_allele_freqs=sapply(output, "[[", "majorFreqs"),
                               minor_allele_freqs=sapply(output, "[[", "minorFreqs")) 
    if (flagRareAlleles) {
        output_table[,"rare_flags"] <- sapply(output, "[[", "minorAlleleFlags")
    }
    return(output_table)
}
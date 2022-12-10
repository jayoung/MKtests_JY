#### excel output


#### addDataAndFormatWorksheet - used in other functions for Excel output
addDataAndFormatWorksheet <- function(myWB, df, sheetName, 
                                      keepNA=TRUE, 
                                      NAcharacter="N.A.", # only used if keepNA==FALSE
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
                                      pop1alias=NULL, pop2alias=NULL,
                                      extraVerbose=FALSE) {
    require(openxlsx)
    ## prep the data frame
    df2 <- df
    # in p-value and NI and alpha columns, replace NA/NaN with "N.A." and infinite values with "Inf". Keep track of which columns I have some characters in, so I can avoid making Excel think they're numeric
    
    columnsContainingCharacter <- list()
    
    for (tempCol in grep("pVal|_NI$|_alpha$",colnames(df2),value=TRUE)) {
        ## I was previously ALWAYS replacing NA/NaN with "N.A." but that messes up Excel formatting when I want to round all the other values in the column. If I switch the NAs to a text value (NAcharacter) I need to do any rounding manually
        if (!keepNA) {
            if (roundSomeColumns) {
                df2[,tempCol] <- round(df2[,tempCol], digits=3)
            }
            df2[which(is.nan(df2[,tempCol])),tempCol] <- NAcharacter
            df2[which(is.na(df2[,tempCol])),tempCol] <- NAcharacter
        }
        df2[which(is.infinite(df2[,tempCol])),tempCol] <- "Inf"
        if(sum(df2[,tempCol] %in% c("Inf", NAcharacter) )>0) {
            columnsContainingCharacter[[tempCol]] <- 1
        }
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
    if(extraVerbose) {cat("here1 colnames are now",colnames(df2),"\n")}
    colnames(df2) <- gsub("_"," ", colnames(df2))
    if(extraVerbose) {cat("here2 colnames are now",colnames(df2),"\n")}
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
    writeData(myWB, sheetName, df2, keepNA=keepNA,
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
        # if we replaced NA values with a text string, Excel gets fussy if we convert to number format, if there were any NAs, but I'm not going to the trouble of testing to see if there were
        if (keepNA) {
            columnsToRound <- grep(" pVal| NI| alpha", colnames(df2))
            if (extraVerbose) { cat("columnsToRound ",columnsToRound,"\n") }
            if (length(columnsToRound)>0) {
                addStyle(myWB, sheetName, createStyle(numFmt="0.000"), stack=TRUE, 
                         cols=columnsToRound, 
                         rows=0:dim(df2)[1]+1, 
                         gridExpand=TRUE)
            }
        }
    }
    
    # some columns need to be converted to text
    if (length(columnsContainingCharacter)>0) {
        for (tempCol in names(columnsContainingCharacter)) {
            cat("converting column",tempCol,"to text\n")
            addStyle(myWB, sheetName, createStyle(numFmt="TEXT"), stack=TRUE, 
                     cols=tempCol, 
                     rows=0:dim(df2)[1]+1, 
                     gridExpand=TRUE)
        }
    }
    
}

### use openxlsx to write an excel file for MKtable and positionTable (two sheets)
writeOneMKtestToExcel <- function(finalOutputTable, positionTable, 
                                  outfileMK,
                                  pop1alias=NULL, pop2alias=NULL,
                                  keepNA = TRUE,
                                  NAcharacter="N.A." # only used if keepNA==FALSE
                                  ) {
    require(openxlsx)
    cat("    saving tables to Excel file\n")
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
    ### the MK results
    addDataAndFormatWorksheet(wb, 
                              finalOutputTable, "MKtable", 
                              rotateColNames=90, 
                              colNamesForVerticalBorder=colsForVerticalBorderResultsTables,
                              pop1alias=pop1alias, pop2alias=pop2alias,
                              headerRowHeight=180, 
                              keepNA=keepNA, NAcharacter=NAcharacter)
    ### the table showing what's going on at each position
    addDataAndFormatWorksheet(wb, 
                              positionTable, "nucleotidePositions", 
                              rotateColNames=90, numColumnsToFreeze=3,
                              colNamesForVerticalBorder=colsForVerticalBorderPositionTables, 
                              rowIndicesForHorizontalBorder=rowIndicesForHorizontalBorder, 
                              rowIndicesColorList=list(green3=polyPositionsToColor, 
                                                       red=fixedPositionsToColor),
                              pop1alias=pop1alias, pop2alias=pop2alias,
                              headerRowHeight=150)
    saveWorkbook(wb, outfileMK, overwrite = TRUE) # , keepNA = TRUE
    rm(wb)
}
    
### combineMKresults: takes MK results from several alignments and makes a single output file / table
combineMKresults <- function(MKresultList, outFile=NULL, outDir=NULL, 
                             keepNA=TRUE, NAcharacter="N.A.",
                             pop1alias=NULL, pop2alias=NULL,
                             getGeneNames=FALSE, geneNameFile="riniTable2_geneOrder.txt", 
                             rowBordersEachGene=FALSE,
                             roundSomeColumns=TRUE,
                             extraVerbose=FALSE) {
    ## once in a while I want to combine but NOT save an excel file
    if(!is.null(outFile)) {
        require(openxlsx)
        #### check outDir is present
        if(!is.null(outDir)) {
            if (!dir.exists(outDir)) { dir.create(outDir) }
            outFile <- paste(outDir,outFile,sep="/")
        }
    }
    #### check they have names
    if(is.null(names(MKresultList))) {
        stop("\n\nERROR - the results list must have names\n\n")
    }
    if(length(unique(names(MKresultList))) != length(MKresultList) ) {
        stop("\n\nERROR - each item in the results list must have a unique name\n\n")
    }
    #cat("keepNA is",keepNA,"\n")
    
    #### if the list has a mix of unpolarized and polarized results, I need to do something more before I can rbind.
    summaries <- lapply(MKresultList, "[[", "summary")
    
    ## first we choose a summary table whose colnames are what we'll aim for in all the other summary tables. which.max arbitrarily chooses the first if there are ties
    maxNumCols <- which.max(sapply(summaries, ncol))
    summaryWithMostColumns <- summaries [[ maxNumCols ]]
    colnamesWeWant <- colnames(summaryWithMostColumns)
    
    ## then we make sure each summary table has those columns, in that order. We fill in missing columns with "". We will get weird Excel warnings (number stored as text) but we'll have to ignore that
    summariesFixed <- lapply(summaries, function(x) {
        # maybe colnames are already identical and we don't need to do anything
        if (identical(colnamesWeWant, colnames(x))) {return(x)}
        
        # maybe the same colnames exist but we need to reorder them
        if (identical(sort(colnamesWeWant), sort(colnames(x)))) {
            x <- x[,colnamesWeWant]
            return(x)
        }
                
        # maybe there are the same NUMBER of columns but they have different names (I'm not testing explicitly for different names, but the two conditions above should have already picked up situations where the names are the same)
        if (length(colnamesWeWant) == length(colnames(x))) {
            stop("\n\nERROR - cannot combine these summary tables, because they have different column names. The code is not (yet?) set up to handle this\n\n")
        }
        if (!keepNA & !is.null(outFile)) {
            stop("\n\nERROR - haven't figured out the correct way to code this situation. We are trying to merge unpolarized and polarized results into a single output file, and we are trying to use the keepNA=FALSE option while exporting to Excel. It's acting weirdly and I haven't done figured out why yet\n")
        }
        
        # maybe some columns are totally missing
        missingColnames <- setdiff(colnamesWeWant, colnames(x))
        if (length(missingColnames)>1) {
            for (missingCol in missingColnames) {
                x[,missingCol] <- NA
                #keepNA <<- FALSE  ## special assignment operator so we affect the global variable for the rest of the function
                ## if we do this, we'll want to mess with the Excel formatting 
                # keepNA=TRUE, NAcharacter="N.A.",
            }
            x <- x[,colnamesWeWant]
            return(x)
        }
    })
    #cat("keepNA is",keepNA,"\n")
    #### combine results into a single table
    results_df <- do.call("rbind", summariesFixed)
    
    if(extraVerbose) { cat("colnames are now", colnames(results_df), "\n") }
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
    
    ### maybe I don't want the Excel file:
    if(is.null(outFile)) {
        return(results_df)
    }
    
    ####### save output to excel file
    results_df_forOutput <- results_df
    wb <- createWorkbook()
    addDataAndFormatWorksheet(wb, keepNA=keepNA, NAcharacter=NAcharacter,
                              results_df_forOutput, "MKtests", 
                              rotateColNames=90, 
                              colNamesForVerticalBorder=colsForVerticalBorderResultsTables, 
                              pop1alias=pop1alias, pop2alias=pop2alias,
                              headerRowHeight=180,
                              roundSomeColumns=roundSomeColumns,
                              extraVerbose=extraVerbose)
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
    #cat("    adding position tables\n")
    positionsTables <- lapply(MKresultList, "[[", "positions")
    names(positionsTables) <- names(MKresultList)
    for (resultName in names(positionsTables)) {
        #cat("    adding position table for",resultName,"\n")
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
        addDataAndFormatWorksheet(wb, keepNA=FALSE,
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


#### plotting 

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
    
    ## check all colnames present
    colNamesNotPresent <- sapply(colNames, function(x) {! x %in% colnames(df) })
    if(sum(colNamesNotPresent)>0) {
        missingColNames <- names(colNamesNotPresent)[which(colNamesNotPresent)]
        cat("\n\nERROR - not all colNames are present in the data frame. The missing ones are:",missingColNames,"\n\n")
        stop()
    }
    
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
    
    ## synon fixed 
    if (sum(df[,colNames[["Ds"]] ]>0, na.rm=TRUE)) {
        Ds_table <- df[,c("codon",colNames[["Dn"]],colNames[["Ds"]])] 
        Ds_table <- Ds_table[which(Ds_table[,colNames[["Ds"]] ]>0),]
        # for y1 we start at pop1_vs_pop2_Dn, so that values will be stacked for codons that have both Dn and Ds
        segments(x0=Ds_table[,"codon"], y0=Ds_table[,colNames[["Dn"]]],
                 x1=Ds_table[,"codon"], y1=Ds_table[,colNames[["Ds"]] ],
                 col=myColors[["Ds"]], lwd=myLwd)
    }
    ## non-synon fixed (plot after synon, to make sure we see them)
    if (sum(df[, colNames[["Dn"]] ]>0, na.rm=TRUE)) {
        Dn_table <- df[,c("codon",colNames[["Dn"]])] 
        Dn_table <- Dn_table[which(Dn_table[,colNames[["Dn"]] ]>0),]
        segments(x0=Dn_table[,"codon"], y0=0,
                 x1=Dn_table[,"codon"], y1=Dn_table[,colNames[["Dn"]] ],
                 col=myColors[["Dn"]], lwd=myLwd)
    }
    ## synon poly
    if (sum(df[, colNames[["Ps"]] ]>0, na.rm=TRUE)) {
        Ps_table <- df[which(df[, colNames[["Ps"]] ]>0),]
        segments(x0=Ps_table[,"codon"], y0=(0-Ps_table[,colNames[["Pn"]] ]),
                 x1=Ps_table[,"codon"], y1=(0-Ps_table[,colNames[["Ps"]] ]),
                 col=myColors[["Ps"]], lwd=myLwd)
    }
    ## non-synon poly
    if (sum(df[,colNames[["Pn"]]]>0, na.rm=TRUE)) {
        Pn_table <- df[which(df[, colNames[["Pn"]] ]>0),]
        segments(x0=Pn_table[,"codon"], y0=0,
                 x1=Pn_table[,"codon"], y1=(0-Pn_table[, colNames[["Pn"]] ]),
                 col=myColors[["Pn"]], lwd=myLwd)
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
                            plotPolarizedPop1=FALSE, plotPolarizedPop2=FALSE, 
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

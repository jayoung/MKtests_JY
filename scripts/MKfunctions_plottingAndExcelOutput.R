#### excel output

### use openxlsx to write an excel file for MKtable and positionTable (two sheets)
writeOneMKtestToExcel <- function(finalOutputTable, positionTable, 
                                  outfileMK,
                                  pop1alias=NULL, pop2alias=NULL) {
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
    addDataAndFormatWorksheet(wb, 
                              finalOutputTable, "MKtable", 
                              rotateColNames=90, 
                              colNamesForVerticalBorder=colsForVerticalBorderResultsTables,
                              pop1alias=pop1alias, pop2alias=pop2alias,
                              headerRowHeight=180)
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

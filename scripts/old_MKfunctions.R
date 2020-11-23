
####### MKoutput works on a table of changes at each nucleodide/codon like positionTable
## output is modelled on Rini's table 2
MKoutput_forMelSim <- function(df) {
    changeTotals <- apply( df[,grep("[DP][ns]", colnames(df), value=TRUE, perl=TRUE)], 2, sum, na.rm=TRUE)
    ## numbers for regular MK test (just mel-sim, no polarization):
    melSimNums <- list(Dn=changeTotals["melSim_Dn"], Ds=changeTotals["melSim_Ds"], Pn=changeTotals["mel_Pn"], Ps=changeTotals["mel_Ps"] )
    ## numbers for polarized MK test (just changes on mel branch versus mel polymorphisms):
    polarizedMelNums <- list(Dn=changeTotals["melPolarized_Dn"], Ds=changeTotals["melPolarized_Ds"], Pn=changeTotals["mel_Pn"], Ps=changeTotals["mel_Ps"] )
    output <- data.frame( melSim_pVal=do.call("doChiSq", melSimNums),
                          melSim_FET_pVal=do.call("doFisher", melSimNums),
                          melSim_Dn=changeTotals["melSim_Dn"],
                          melSim_Ds=changeTotals["melSim_Ds"],
                          melSim_Pn=changeTotals["mel_Pn"],
                          melSim_Ps=changeTotals["mel_Ps"],
                          melSim_NI=do.call("neutralityIndex",melSimNums),
                          polarizedMel_pVal=do.call("doChiSq", polarizedMelNums),
                          polarizedMel_FET_pVal=do.call("doFisher", polarizedMelNums),
                          polarizedMel_Dn=changeTotals["melPolarized_Dn"], 
                          polarizedMel_Ds=changeTotals["melPolarized_Ds"], 
                          polarizedMel_NI=do.call("neutralityIndex",polarizedMelNums), 
                          polarizedSim_Dn=changeTotals["simPolarized_Dn"], 
                          polarizedSim_Ds=changeTotals["simPolarized_Ds"],
                          total_polarized_Dn=changeTotals["total_polarized_Dn"], 
                          total_polarized_Ds=changeTotals["total_polarized_Ds"])
    return(output)
}




saveExcelForCombinedMKoutput <- function(MKresult_df, 
                                         outFile, 
                                         rowBordersEachGene=FALSE, 
                                         pop1alias=NULL, pop2alias=NULL) {
    require(openxlsx)
    ####### save output to excel file
    results_df_forOutput <- MKresult_df
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
        rm(tempRowsToAddBorder)
    }
    saveWorkbook(wb, outFile, overwrite = TRUE)
}


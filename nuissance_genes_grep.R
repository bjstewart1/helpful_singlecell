#core exclude genes
core_exclude_grep <- function(charvect){
hkGeneREGEX='^(DNAJ[ABC]|EIF[0-9]|RPL[0-9]|RPS[0-9]|RPN1|POLR[0-9]|SNX[0-9]|HSP[AB][0-9]|H1FX|H2AF[VXYZ]|PRKA|NDUF[ABCSV]|ATP[0-9]|PSM[ABCDEFG][0-9]|UBA[0-9]|UBE[0-9]|USP[0-9]|TXN)'
coreExcludeGenes = unique(c(grep('\\.[0-9]+',charvect,value=TRUE), #Poorly characterised
                            grep('MALAT1',charvect,value=TRUE), #Contamination or highly expressed poorly characterised
                            grep('^MT-',charvect,value=TRUE), #Mitochondria
                            grep("XIST", charvect, value = TRUE), #F gender
                            grep(hkGeneREGEX,charvect,value=TRUE) #Housekeeping genes
))
core_exclude <- !charvect %in% coreExcludeGenes
return(core_exclude)
}

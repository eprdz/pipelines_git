## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- echo=FALSE--------------------------------------------------------------
doAll=FALSE

## ---- eval=doAll, results="hide"----------------------------------------------
#  pacman::p_load(edgeR, enriquepresa, ReportingTools, SummarizedExperiment)
#  wd = getwd()

## ---- eval=doAll--------------------------------------------------------------
#  data(PRJNA517180, package = "enriquepresa")
#  grupos = colData(PRJNA517180)[,"type"]
#  
#  x = DGEList(counts = assay(PRJNA517180), group = grupos)

## ---- eval=doAll--------------------------------------------------------------
#  keep = rowSums(cpm(x) > 1) >= 2
#  x <- x[keep, , keep.lib.sizes = FALSE]

## ---- eval=doAll--------------------------------------------------------------
#  dge.c = estimateCommonDisp(x)
#  dge.t = estimateTagwiseDisp(dge.c)

## ---- eval=doAll--------------------------------------------------------------
#  et.t = exactTest(dge.t)
#  df = topTags(et.t, n=nrow(et.t$table))
#  df = df$table

## ---- eval=doAll--------------------------------------------------------------
#  SYMBOL = rowData(PRJNA517180)[rownames(df),"SYMBOL"]
#  ENTREZ_URLs = sapply(rowData(PRJNA517180)[rownames(df),"ENTREZID"], entrezurl)
#  ENSEMBL_URLs = sapply(rowData(PRJNA517180)[rownames(df),"ENSEMBL"], ensemblurl)
#  df = cbind(SYMBOL, df, ENTREZ_URLs, ENSEMBL_URLs)

## ---- eval=doAll, warning=FALSE, message=FALSE--------------------------------
#  foutput = "PRJNA517180"
#  htmlRep1 = HTMLReport(shortName = foutput,title = foutput,
#                        reportDirectory = "./reports/differential_expression")
#  
#  publish(df, htmlRep1)
#  finish(htmlRep1)

## ---- eval=doAll--------------------------------------------------------------
#  
#  edgeR_df = df
#  dirReports = paste0(wd,"/reports/")
#  save(edgeR_df, file=paste0(dirReports,"edgeR_df.rda"))

## ---- eval=doAll--------------------------------------------------------------
#  plot_volcano(edgeR_df$FDR, edgeR_df$logFC, edgeR_df$SYMBOL, 0.05, 1)

## ---- eval=doAll--------------------------------------------------------------
#  sum(edgeR_df$FDR < 0.05 & edgeR_df$logFC > 1) #Número de genes sobreexpresados
#  sum(edgeR_df$FDR < 0.05 & edgeR_df$logFC < -1) #Número de genes subexpresados


## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- echo=FALSE--------------------------------------------------------------
doAll=FALSE

## ---- eval=doAll, results="hide"----------------------------------------------
#  pacman::p_load(BiocGenerics, EnrichmentBrowser, enriquepresa, ReportingTools, SummarizedExperiment)

## ---- eval=doAll, warning=FALSE, message=FALSE--------------------------------
#  data(PRJNA517180, package = "enriquepresa")
#  data(edgeR_df, package = "enriquepresa")
#  
#  cel_go_bp = getGenesets(org = "cel", db = "go", go.onto = "BP")
#  cel_go_cc = getGenesets(org = "cel", db = "go", go.onto = "CC")
#  cel_go_mf = getGenesets(org = "cel", db = "go", go.onto = "MF")

## ---- eval=doAll--------------------------------------------------------------
#  PRJNA517180 = PRJNA517180[which(is.na(rowData(PRJNA517180)$ENTREZID)==FALSE),]
#  edgeR_df = edgeR_df[which(is.na(edgeR_df$ENTREZ_URLs)==FALSE),]
#  PRJNA517180 = PRJNA517180[rownames(edgeR_df),]

## ---- eval=doAll--------------------------------------------------------------
#  PRJNA517180 = PRJNA517180[match(unique(rowData(PRJNA517180)$ENTREZID),rowData(PRJNA517180)$ENTREZID),]
#  edgeR_df = edgeR_df[match(unique(edgeR_df$ENTREZ_URLs), edgeR_df$ENTREZ_URLs),]

## ---- eval=doAll--------------------------------------------------------------
#  rownames(PRJNA517180) = rowData(PRJNA517180)[rownames(PRJNA517180),"ENTREZID"]
#  rowData(PRJNA517180) = edgeR_df

## ---- eval=doAll--------------------------------------------------------------
#  colnames(rowData(PRJNA517180)) = c("SYBMOL", "FC", "logCPM", "PValue", "ADJ.PVAL", "ENTREZ_URLs", "ENSEMBL_URLs")
#  colnames(colData(PRJNA517180)) = c("Run", "GROUP")

## ---- eval=doAll--------------------------------------------------------------
#  PRJNA517180.gobp = sbea(method = "ora", se = PRJNA517180, gs = cel_go_bp, perm = 0, alpha = 0.05)
#  PRJNA517180.gocc = sbea(method = "ora", se = PRJNA517180, gs = cel_go_cc, perm = 0, alpha = 0.05)
#  PRJNA517180.gomf = sbea(method = "ora", se = PRJNA517180, gs = cel_go_mf, perm = 0, alpha = 0.05)

## ---- eval=doAll--------------------------------------------------------------
#  PRJNA517180.gobp = gsRanking(PRJNA517180.gobp)
#  PRJNA517180.gocc = gsRanking(PRJNA517180.gocc)
#  PRJNA517180.gomf = gsRanking(PRJNA517180.gomf)
#  rownames(PRJNA517180.gobp) = PRJNA517180.gobp[,"GENE.SET"]
#  rownames(PRJNA517180.gocc) = PRJNA517180.gocc[,"GENE.SET"]
#  rownames(PRJNA517180.gomf) = PRJNA517180.gomf[,"GENE.SET"]

## ---- eval=doAll, warning=FALSE, message=FALSE--------------------------------
#  PRJNA517180_ora.gobp = enriquepresa::dfGO(PRJNA517180.gobp)
#  PRJNA517180_ora.gocc = enriquepresa::dfGO(PRJNA517180.gocc)
#  PRJNA517180_ora.gomf = enriquepresa::dfGO(PRJNA517180.gomf)
#  head(PRJNA517180_ora.gobp) #Un ejemplo de como queda con las URLs
#  
#  foutput = "PRJNA517180_ora.gobp"
#  htmlRep1 = HTMLReport(shortName = foutput,title = foutput,
#                        reportDirectory = "./reports/ora")
#  
#  publish(PRJNA517180_ora.gobp, htmlRep1)
#  finish(htmlRep1)
#  
#  foutput = "PRJNA517180_ora.gocc"
#  htmlRep1 = HTMLReport(shortName = foutput,title = foutput,
#                        reportDirectory = "./reports/ora")
#  
#  publish(PRJNA517180_ora.gocc, htmlRep1)
#  finish(htmlRep1)
#  
#  foutput = "PRJNA517180_ora.gomf"
#  htmlRep1 = HTMLReport(shortName = foutput,title = foutput,
#                        reportDirectory = "./reports/ora")
#  
#  publish(PRJNA517180_ora.gomf, htmlRep1)
#  finish(htmlRep1)

## ---- eval=doAll, warning=FALSE, message=FALSE--------------------------------
#  data(Eset50467, package = "enriquepresa")
#  
#  hsa_go_bp = getGenesets(org = "hsa", db = "go", go.onto = "BP")
#  hsa_go_cc = getGenesets(org = "hsa", db = "go", go.onto = "CC")
#  hsa_go_mf = getGenesets(org = "hsa", db = "go", go.onto = "MF")

## ---- eval=doAll, warning=FALSE, message=FALSE--------------------------------
#  se50467 = makeSummarizedExperimentFromExpressionSet(Eset50467)
#  se50467 = probe2gene(se50467)
#  GROUP = c(rep(0, 6), rep(1, 6))
#  colData(se50467) = cbind(colData(se50467), GROUP)
#  se50467 = deAna(expr = se50467, de.method = "limma")

## ---- eval=doAll--------------------------------------------------------------
#  se50467.gobp = sbea(method = "ora", se = se50467, gs = hsa_go_bp, perm = 0, alpha = 0.05)
#  se50467.gocc = sbea(method = "ora", se = se50467, gs = hsa_go_cc, perm = 0, alpha = 0.05)
#  se50467.gomf = sbea(method = "ora", se = se50467, gs = hsa_go_mf, perm = 0, alpha = 0.05)

## ---- eval=doAll--------------------------------------------------------------
#  se50467.gobp = gsRanking(se50467.gobp)
#  se50467.gocc = gsRanking(se50467.gocc)
#  se50467.gomf = gsRanking(se50467.gomf)
#  rownames(se50467.gobp) = se50467.gobp[,"GENE.SET"]
#  rownames(se50467.gocc) = se50467.gocc[,"GENE.SET"]
#  rownames(se50467.gomf) = se50467.gomf[,"GENE.SET"]

## ---- eval=doAll, warning=FALSE, message=FALSE--------------------------------
#  se50467_ora.gobp = enriquepresa::dfGO(se50467.gobp)
#  se50467_ora.gocc = enriquepresa::dfGO(se50467.gocc)
#  se50467_ora.gomf = enriquepresa::dfGO(se50467.gomf)
#  head(se50467_ora.gobp)
#  
#  foutput = "se50467_ora.gobp"
#  htmlRep1 = HTMLReport(shortName = foutput,title = foutput,
#                        reportDirectory = "./reports/ora")
#  
#  publish(se50467_ora.gobp, htmlRep1)
#  finish(htmlRep1)
#  
#  foutput = "se50467_ora.gocc"
#  htmlRep1 = HTMLReport(shortName = foutput,title = foutput,
#                        reportDirectory = "./reports/ora")
#  
#  publish(se50467_ora.gocc, htmlRep1)
#  finish(htmlRep1)
#  
#  foutput = "se50467_ora.gomf"
#  htmlRep1 = HTMLReport(shortName = foutput,title = foutput,
#                        reportDirectory = "./reports/ora")
#  
#  publish(se50467_ora.gomf, htmlRep1)
#  finish(htmlRep1)


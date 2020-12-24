## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- echo=FALSE--------------------------------------------------------------
doAll=FALSE

## ---- eval=doAll, results="hide"----------------------------------------------
#  pacman::p_load(Biobase, enriquepresa, genefilter, limma, ReportingTools)

## ---- eval=doAll, results="hide"----------------------------------------------
#  data(Eset50467, package="enriquepresa")

## ---- eval=doAll--------------------------------------------------------------
#  grupos = pData(Eset50467)[,"FactorValue..compound."]
#  genefilter_df = rowttests(Eset50467, grupos)
#  colnames(genefilter_df)[2] = "LogFC"
#  head(genefilter_df)

## ---- eval=doAll--------------------------------------------------------------
#  BH = p.adjust(genefilter_df$p.value, method = "BH")
#  Bonferroni = p.adjust(genefilter_df$p.value, method = "bonferroni")
#  genefilter_df = cbind(genefilter_df, BH, Bonferroni)
#  genefilter_df = genefilter_df[order(BH),] #Se ordenan las filas en orden decreciente según el adj. p-value según el método BH

## ---- eval=doAll--------------------------------------------------------------
#  fd = fData(Eset50467)
#  entrezids = sapply(rownames(genefilter_df), function(x){fd[fd$PROBEID==x,"ENTREZID"]})
#  ENTREZ_URL = sapply(entrezids, entrezurl)
#  ensemblids = sapply(rownames(genefilter_df), function(x){fd[fd$PROBEID==x,"ENSEMBL"]})
#  ENSEMBL_URL = sapply(ensemblids, ensemblurl)
#  SYMBOL = sapply(rownames(genefilter_df), function(x){fd[fd$PROBEID==x,"SYMBOL"]})
#  
#  genefilter_df = cbind(SYMBOL, genefilter_df, ENTREZ_URL, ENSEMBL_URL)
#  head(genefilter_df)

## ---- eval=doAll--------------------------------------------------------------
#  foutput = "GSE50467_rowttests"
#  htmlRep1 = HTMLReport(shortName = foutput,title = foutput,
#                        reportDirectory = "./reports/differential_expression")
#  
#  publish(genefilter_df, htmlRep1)
#  finish(htmlRep1)

## ---- eval=doAll--------------------------------------------------------------
#  save(genefilter_df, file="./reports/genefilter_df.rda", compress = "xz")

## ---- eval=doAll--------------------------------------------------------------
#  grupos= pData(Eset50467)[,"FactorValue..compound."]
#  design = model.matrix(~0+grupos)
#  colnames(design) = c("doxycycline", "control")
#  design

## ---- eval=doAll--------------------------------------------------------------
#  fit = lmFit(Eset50467, design)
#  cont.matrix <- makeContrasts(doxycycline - control, levels=design)

## ---- eval=doAll--------------------------------------------------------------
#  fit2 <- contrasts.fit(fit, cont.matrix)
#  fit2 <- eBayes(fit2)
#  tT <- topTable(fit2, coef = 1, adjust = "BH", number = nrow(Eset50467))
#  tT2 <- topTable(fit2, coef = 1, adjust = "bonferroni", number = nrow(Eset50467))

## ---- eval=doAll--------------------------------------------------------------
#  entrezurls = sapply(tT$ENTREZID, entrezurl)
#  ensemblurls = sapply(tT$ENSEMBL, ensemblurl)
#  limma_df = data.frame(SYMBOL = tT$SYMBOL, statistic = tT$t, LogFC=tT$logFC, p.value=tT$P.Value,
#                  BH=tT$adj.P.Val, Bonferroni=tT2$adj.P.Val, ENTREZ_URL=entrezurls, ENSEMBL_URL=ensemblurls)
#  rownames(limma_df) = tT$PROBEID
#  head(limma_df)

## ---- eval=doAll--------------------------------------------------------------
#  foutput = "GSE50467_limma"
#  htmlRep1 = HTMLReport(shortName = foutput,title = foutput,
#                        reportDirectory = "./reports/differential_expression")
#  
#  publish(limma_df, htmlRep1)
#  finish(htmlRep1)

## ---- eval=doAll--------------------------------------------------------------
#  save(limma_df, file="./reports/limma_df.rda", compress = "xz")

## ---- eval=doAll--------------------------------------------------------------
#  plot_volcano(genefilter_df$BH, genefilter_df$LogFC, genefilter_df$SYMBOL, min_p = 0.05, diff = 1)

## ---- eval=doAll--------------------------------------------------------------
#  2^(as.numeric(limma_df[limma_df$SYMBOL=="VEGFA","LogFC"]))
#  1/(2^(as.numeric(limma_df[limma_df$SYMBOL=="EGFR","LogFC"])))


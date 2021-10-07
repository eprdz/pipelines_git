## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- echo=FALSE--------------------------------------------------------------
doAll=FALSE

## ---- eval=doAll, results="hide"----------------------------------------------
#  pacman::p_load(Biobase, GEOquery, limma, tidyr, AnnotationDbi, BiocGenerics, hgug4112a.db)

## ---- eval=doAll, results="hide"----------------------------------------------
#  wd = getwd()
#  system("mkdir dirData")
#  dirData = paste0(wd,"/dirData/")
#  setwd(dirData)
#  GEOquery::getGEOSuppFiles("GSE50467")
#  setwd("GSE50467")
#  system("tar xvf GSE50467_RAW.tar")
#  system("gzip -d *.gz")
#  x=limma::read.maimages(dir(".","txt"),"agilent",green.only=TRUE,other.columns="gIsWellAboveBG")
#  GSE50467raw = x
#  save(GSE50467raw, file=paste0(dirData,"GSE50467raw.rda"))
#  setwd(wd)

## ---- eval=doAll, message=FALSE-----------------------------------------------
#  load(paste0(dirData,"GSE50467raw.rda"))
#  annot = AnnotationDbi::select(hgug4112a.db, keys=GSE50467raw$genes[,"ProbeName"],
#  column = c("ENTREZID","ENSEMBL"), keytype="PROBEID")
#  annot = annot[!is.na(annot[,"ENTREZID"]),]
#  uniq_probe = match(unique(annot[,1]),annot[,1])
#  annot1 = annot[uniq_probe,]
#  uniq_entrez = match(unique(annot1[,2]), annot1[,2])
#  annot2 = annot1[uniq_entrez,]
#  annot2 = na.omit(annot2)
#  GSE50467raw = GSE50467raw[match(annot2[,1],GSE50467raw$genes[,"ProbeName"]),]
#  rownames(GSE50467raw) = annot2[,1]
#  save(GSE50467raw, file=paste0(dirData,"GSE50467raw.rda"))

## ---- warning = FALSE, eval=doAll---------------------------------------------
#  load(paste0(dirData,"GSE50467raw.rda"))
#  n.arrays = c(1:ncol(GSE50467raw))
#  invisible(lapply(n.arrays,function(x) {plotMA(GSE50467raw, array=x)}))
#  plotDensities(GSE50467raw,log=TRUE,legend=FALSE, main="Histograma de las muestras sin normalizar")
#  boxplot(GSE50467raw$E, xlab="Muestras", main="Boxplot de las muestras sin normalizar", xaxt="n")

## ---- eval=doAll, results="hide"----------------------------------------------
#  GSE50467=backgroundCorrect(GSE50467raw, method="normexp")
#  GSE50467=normalizeBetweenArrays(GSE50467, method="quantile")
#  save(GSE50467, file=paste0(dirData,"GSE50467.rda"))

## ---- eval=doAll--------------------------------------------------------------
#  load(paste0(dirData,"GSE50467.rda"))
#  setwd(dirData)
#  system("wget https://www.ebi.ac.uk/arrayexpress/files/E-GEOD-50467/E-GEOD-50467.sdrf.txt")
#  fenodata = read.csv("E-GEOD-50467.sdrf.txt",sep="\t",header=TRUE)
#  setwd(wd)
#  fenodata = fenodata[order(fenodata$Source.Name),]
#  pd = new("AnnotatedDataFrame", data = fenodata)
#  rownames(pd) = fenodata$Source.Name

## ----echo=FALSE, eval=doAll---------------------------------------------------
#  rownames(pd)

## ---- eval=doAll--------------------------------------------------------------
#  cols = c()
#  colnames(GSE50467$E) = as.vector(lapply(colnames(GSE50467$E), function(x){c(cols, paste0(substr(x, 1, 10), " 1"))}))

## ---- eval=doAll--------------------------------------------------------------
#  setwd(dirData)
#  system("wget https://www.ebi.ac.uk/arrayexpress/files/E-GEOD-50467/E-GEOD-50467.idf.txt")
#  experimentdata = read.csv("E-GEOD-50467.idf.txt",header = FALSE, sep="\t")
#  setwd(wd)
#  experimentdata2 = t(tidyr::unite(experimentdata, "data", 2:8, sep=" "))
#  exp.names = experimentdata2[1,]
#  exp.list = as.list(experimentdata2[-1,])
#  names(exp.list) = exp.names
#  MIAME = MIAME(name=exp.list$`Publication Author List`, lab = exp.list$`Person Address`, contact = exp.list$`Person Email`, title = exp.list$`Investigation Title`, abstract=exp.list$`Experiment Description`, url = paste0("http://dx.doi.org/",substr(exp.list$`Publication DOI`, 1,22)), pubMedIds = substr(exp.list$`Pubmed ID`, 1, 8), other = list(ExtraInfo = 'MIAME created from list with experimental data.'))

## ---- eval=doAll--------------------------------------------------------------
#  MIAME

## ---- eval=doAll--------------------------------------------------------------
#  rownames(annot2) = annot2[,1]
#  fD = new("AnnotatedDataFrame", data = annot2)

## ---- eval=doAll--------------------------------------------------------------
#  Exp.set = new("ExpressionSet",exprs=GSE50467$E,phenoData=pd,experimentData = MIAME, featureData=fD, annotation = "hgug4112a.db")
#  save(Exp.set, file=paste0(wd,"/Eset50467.rda"))

## ---- warning=FALSE, eval=doAll-----------------------------------------------
#  n.arrays = c(1:ncol(GSE50467))
#  invisible(lapply(n.arrays,function(x) {plotMA.EList(GSE50467, array=x)}))
#  plotDensities(GSE50467,log=TRUE,legend=FALSE, main="Histograma de las muestras normalizadas")
#  boxplot(GSE50467$E, xlab="Muestras", main="Boxplot de las muestras normalizadas", xaxt="n")

## ---- eval=doAll, message=FALSE-----------------------------------------------
#  final_annot = AnnotationDbi::select(hgug4112a.db, keys=fData(Exp.set)$PROBEID,
#                                column=c("ENTREZID","ENSEMBL", "SYMBOL"), keytype="PROBEID")
#  fData(Exp.set) = final_annot[match(featureNames(Exp.set),final_annot$PROBEID),]
#  head(fData(Exp.set))
#  save(Exp.set, file=paste0(wd,"/Eset50467.rda"))


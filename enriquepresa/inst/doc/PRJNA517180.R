## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- echo=FALSE--------------------------------------------------------------
doAll=FALSE

## ---- eval=doAll--------------------------------------------------------------
#  wd = getwd()

## ---- eval=doAll--------------------------------------------------------------
#  pacman::p_load(Rsamtools, GenomicFeatures, GenomicAlignments)
#  setwd("~/ncbi/public")
#  gtfFile = "Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/genes.gtf"
#  txdb = makeTxDbFromGFF(gtfFile, format="gtf")
#  genes = exonsBy(txdb, by="gene")
#  
#  setwd("~/ncbi/public/aligned")
#  dirActualData =  paste(getwd(),"/",sep="")
#  sampleTableSingle = read.table("bamfiles.txt")
#  fls = paste(dirActualData,sampleTableSingle[,1],sep="")
#  bamLst = BamFileList(fls, index=character(),yieldSize=100000,obeyQname=TRUE)
#  PRJNA517180 = summarizeOverlaps(features = genes,read=bamLst,
#     mode="Union",
#        singleEnd=FALSE,
#        ignore.strand=TRUE,
#        fragments=TRUE)
#  Run = c("SRR8489772", "SRR8489773", "SRR8489774", "SRR8489775", "SRR8489776", "SRR8489777", "SRR8489778", "SRR8489779")
#  
#  type = c(rep(0, 4), rep(1, 4))
#  type = factor(type, levels=0:1, labels=c("L4440", "ztf_11_KD"))
#  colData(PRJNA517180) = DataFrame(Run, type)

## ---- eval=doAll--------------------------------------------------------------
#  pacman::p_load(AnnotationDbi, org.Ce.eg.db)
#  annot = AnnotationDbi::select(org.Ce.eg.db, keys=rownames(assay(PRJNA517180)),
#                                column = c("SYMBOL", "ENTREZID", "ENSEMBL"), keytype="WORMBASE")
#  uniq = match(rownames(assay(PRJNA517180)),annot[,1])
#  rowData(PRJNA517180)=annot[uniq,]
#  save(PRJNA517180, file=paste0(wd,"/PRJNA517180.rda"))


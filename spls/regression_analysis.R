#!/usr/bin/env Rscript

####################################
####     Regression analysis    ####
####################################

pacman::p_load(mixOmics, ropls)

parameters <- commandArgs(trailingOnly=TRUE)

counts.file <- as.character(parameters[1]) # Omic dataset
meta.file <- as.character(parameters[2]) # Metadata file
y <- as.character(parameters[3]) # Continuous variables to study (logFCcfu)
descr.column <- as.character(parameters[4]) # What is going to be analysed - Genera, metabolites or IDs from KEGG
ncomp <- as.numeric(as.character(parameters[5])) # Number of components to use (2)
omic <- as.character(parameters[6]) # 16S rRNA analysis, metabolomics, proteomics or metagenomics
model <- as.character(parameters[7]) # dynamic model or predictive model

##########
## Data selection and arrange
dat <- read.table(file=counts.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE)
DESCRIPTION <- dat[, descr.column]; names(DESCRIPTION) <- rownames(dat) <- dat[, 1]

meta <- read.table(file=meta.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE)
rownames(meta) <- meta[, 1]
mysamples <- intersect(rownames(meta), colnames(dat))
dat <- as.matrix(dat[, mysamples])
meta <- meta[mysamples, ]

Y <- meta[,y]
X <- t(dat)

##########
## sPLS-reg

# Grid of possible numbers of variables that will be tested for each component
if (omic == "metagenomics" | omic == "proteomics") {
  list.keepX <- c(seq(10, 1010, 50))
} else {
  list.keepX <- c(seq(10, nrow(dat), 5))
  
}

tune.spls <- tune.spls(X, Y, ncomp=ncomp, validation='Mfold', folds=5, 
                           progressBar=TRUE,
                           test.keepX=list.keepX, nrepeat=50)

# The optimal number of features to select (per component):
#tune.spls$choice.keepX

# The optimal number of components
#tune.spls$choice.ncomp$ncomp

# We include these parameters in our final sPLS-reg model:

choice.ncomp <- ncomp
choice.keepX <- tune.spls$choice.keepX[1:choice.ncomp]

# Applying sPLS-reg
spls.res <- mixOmics::spls(X, Y, ncomp=choice.ncomp, keepX=choice.keepX, mode = "regression")
save(spls.res, file="spls.res.RDa")

# Assessing performance of sPLS-reg
perf.spls <- perf(spls.res, validation="Mfold", folds=5, progressBar=TRUE, auc=TRUE, nrepeat=50, criterion="all")
pdf(file="spls.ncomp.pdf")
plot(perf.spls, col=color.mixo(1:3), sd=TRUE, legend.position="horizontal")
dev.off()

# Final selection of features can be output, along with their weight coefficient
# (most important based on their aboslute value) and their frequency in models:
variables1 = c()
variables2 = c()
for(ncomp in 1:choice.ncomp){
  ind.match = match(selectVar(spls.res, comp = ncomp)$X$name, 
                    names(perf.spls$features$stable.X[[ncomp]]))
  Freq = as.numeric(perf.spls$features$stable.X[[ncomp]][ind.match])
  vars.comp = data.frame(selectVar(spls.res, comp = ncomp)$X$value, Freq)
  vars.comp = cbind(ID = rownames(vars.comp), vars.comp)
  if (ncomp == 1) variables1 = as.character(vars.comp[vars.comp$Freq >= 0, 1])
  if (ncomp == 2) variables2 = as.character(vars.comp[vars.comp$Freq >= 0, 1])
  write.table(vars.comp, file=paste("vars.comp", ncomp, "tsv", sep="."), sep="\t", row.names=FALSE, quote=FALSE)
}

# Sample Plots
levels(spls.res$Y) = c("MRE decrease", "MRE increase")
png(paste0("spls_",omic,model,".png"))
plotIndiv(spls.res, comp=c(1,2), pch = c(16,16),group=spls.res$Y, ind.names=F, legend=F,
          title=paste0('sPLS-DA - ', omic, " - ", model,")"), star = T,
          X.label = paste0("Comp 1: ", round(spls.res$explained_variance$X[1], digits = 2)*100,"% Expl. variance"),
          Y.label = paste0("Comp 2: ", round(spls.res$explained_variance$X[2], digits = 2)*100,"% Expl. variance"),
          legend.position = "bottom")
dev.off()

###########
## Assessment of a PLS model with selected variables by sPLS-reg using ropls
if (omic != "proteomics" & omic != "metagenomics" & omic != "16S rRNA analysis") {
  variables = unique(c(variables1, variables2))} else variables = variables1
X = X[,colnames(X) %in% unique(variables)]
storage.mode(X) = "numeric"
pls = ropls::opls(X, Y, predI = 1)

#  VIP values
vips = ropls::getVipVn(pls)
length(vips)
X = X[,colnames(X) %in% names(which(vips>1))]
X2 = t(X)
X2 = cbind(variable = rownames(X2), X2)
write.table(X2, file=paste("do_regression.by.vips", gsub(x = omic, pattern = " ", "."),
                           gsub(x = model, pattern = " ", "."), "tsv", sep="."), sep="\t", quote = F, row.names = F)

# Quality metrics of the model (R2, Q2)
a = as.data.frame(ropls::getSummaryDF(pls))
a = a[, -c(ncol(a)-1, ncol(a))]
write.table(a, file="info_regression.tsv", sep="\t", row.names = F, quote = F)

# Spearman correlation test and enrichment done separately
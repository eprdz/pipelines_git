#!/usr/bin/env Rscript

#######################################
####     Discriminant analysis     ####
#######################################

pacman::p_load(mixOmics, ggplot2, gplots, DiscriMiner, clusterProfiler)

parameters <- commandArgs(trailingOnly=TRUE)

counts.file <- as.character(parameters[1]) # Omic dataset
meta.file <- as.character(parameters[2]) # Metadata file
myfactors <- as.character(parameters[3]) # Factor of metadata to study (Groups)
mylevels <- as.character(parameters[4]) # Levels of the factor to study (Decrease, increase)
descr.column <- as.character(parameters[5]) # What is going to be analysed - Genera, metabolites or IDs from KEGG
ncomp <- as.numeric(as.character(parameters[6])) # Number of components to use (2)
omic <- as.character(parameters[7]) # 16S rRNA, metabolomics, metagenomics or proteomics
model <- as.character(parameters[8]) # dynamic model or predictive model

##########
## Data selection and arrange
dat <- read.table(file=counts.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE)
DESCRIPTION <- dat[, descr.column]; names(DESCRIPTION) <- rownames(dat) <- dat[, 1]

meta <- read.table(file=meta.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE)
rownames(meta) <- meta[, 1]
mysamples <- intersect(rownames(meta), colnames(dat))
dat <- as.matrix(dat[, mysamples])
meta <- meta[mysamples, ]

myfactors.l <- strsplit(myfactors, split=",")
myfactor1 <- myfactors.l[[1]][1]
if(length(myfactors.l[[1]])==2) myfactor2 <- myfactors.l[[1]][2] else myfactor2 <- ""

if(myfactor2!="") GROUPS <- paste(meta[, myfactor1], meta[, myfactor2], sep=".") else GROUPS <- as.character(meta[, myfactor1])
names(GROUPS) <- mysamples

mylevels <- unlist(strsplit(mylevels, split=","))
mylevels <- mylevels[which(is.element(mylevels, unique(GROUPS)))]
ind <- which(is.element(GROUPS, mylevels))

dat <- dat[, ind]
GROUPS <- GROUPS[ind]
meta <- meta[ind, ]

Y <- as.factor(GROUPS)
X <- t(dat)

##########
## sPLS-DA
# Grid of possible numbers of variables that will be tested for each component
if (omic == "metagenomics" | omic == "proteomics") {
  list.keepX <- c(seq(10, 1010, 50))
} else {
  list.keepX <- c(seq(10, ncol(X), 5))
  
}

# Testing the error of different number of variables
tune.splsda <- tune.splsda(X, Y, ncomp=ncomp, validation='Mfold', folds=5, 
                           progressBar=TRUE,
                           test.keepX=list.keepX, nrepeat=50)

# The optimal number of features to select (per component):
#tune.splsda$choice.keepX

# The optimal number of components
#tune.splsda$choice.ncomp$ncomp

# We include these parameters in our final sPLS-DA model:

choice.ncomp <- ncomp
choice.keepX <- tune.splsda$choice.keepX[1:choice.ncomp]


# Applying sPLS-DA
splsda.res <- mixOmics::splsda(X, Y, ncomp=choice.ncomp, keepX=choice.keepX)
save(splsda.res, file="splsda.res.RDa")


# Assessing performance of sPLS-DA
perf.splsda <- perf(splsda.res, validation="Mfold", folds=5, progressBar=TRUE, auc=TRUE, nrepeat=50)
pdf(file="splsda.ncomp.pdf")
plot(perf.splsda, col=color.mixo(1:3), sd=TRUE, legend.position="horizontal")
dev.off()


# Final selection of features can be output, along with their weight coefficient
#(most important based on their aboslute value) and their frequency in models:
variables1 = c()
variables2 = c()
for(ncomp in 1:choice.ncomp){
  ind.match = match(selectVar(splsda.res, comp = ncomp)$name, 
                    names(perf.splsda$features$stable[[ncomp]]))
  Freq = as.numeric(perf.splsda$features$stable[[ncomp]][ind.match])
  vars.comp = data.frame(selectVar(splsda.res, comp = ncomp)$value, Freq)
  vars.comp = cbind(ID = rownames(vars.comp), vars.comp)
  if (ncomp == 1) variables1 = as.character(vars.comp[vars.comp$Freq >= 0, 1])
  if (ncomp == 2) variables2 = as.character(vars.comp[vars.comp$Freq >= 0, 1])
  write.table(vars.comp, file=paste("vars.comp", ncomp, "tsv", sep="."), sep="\t", row.names=FALSE, quote=FALSE)
}


# Sample Plots (PCA and sPLS-DA)
levels(splsda.res$Y) = c("MRE decrease", "MRE increase")
pca = pca(X, center = F, scale = F)
png(paste0("pca_",omic,model,".png"))
plotIndiv(pca, group = splsda.res$Y, ind.names = F, pch = c(16,16),
          title = paste('PCA -', omic, "-", model), star = T, legend = T, legend.position = "bottom")
dev.off()

back = background.predict(splsda.res, comp.predicted = 2) # To color the background according to the group belonging
png(paste0("splsda_",omic,model,".png"))
plotIndiv(splsda.res, comp=c(1,2), rep.space='X-variate', pch = c(16,16),group=splsda.res$Y, ind.names=F, legend=TRUE,
          col.per.group = c("red2", "forestgreen"), title=paste0('sPLS-DA - ', omic, " - ", model), star = T,
          X.label = paste0("Comp 1: ", round(splsda.res$explained_variance$Y[1], digits = 3)*100,"% Expl. variance"),
          Y.label = paste0("Comp 2: ", round(splsda.res$explained_variance$Y[2], digits = 3)*100,"% Expl. variance"),
          background = back, legend.position = "bottom")
dev.off()

###########
## Assessment of a PLS-DA model with selected variables by sPLS-DA using DiscriMiner
if (omic != "proteomics" & omic != "metagenomics" & omic != "16S rRNA analysis") {
  variables = c(variables1, variables2)} else variables = variables1
X = X[,colnames(X) %in% unique(variables)]
storage.mode(X) = "numeric"
discriminer_pls = DiscriMiner::plsDA(variables = X, group = Y, autosel = F, comps = 2, cv = "LOO")

# VIP values
vips = as.data.frame(discriminer_pls$VIP)
vips = cbind(feature = rownames(vips), vips)
vips = vips[order(vips$`Component 1`, decreasing = T),]
write.table(vips, file=paste(omic,"VIPs", model, "tsv", sep="."), sep="\t", quote = F, row.names = F)

# Quality metrics of the model (R2, Q2, error rate)
output = as.data.frame(cbind(discriminer_pls$R2, discriminer_pls$Q2, error_rate = rep(discriminer_pls$error_rate,2)))
output = cbind(comps = rownames(output), output)
write.table(output, file=paste(omic,"r2q2error", model, "tsv", sep="."), sep="\t", quote = F, row.names = F)
vips = discriminer_pls$VIP
X = X[,colnames(X) %in% names(which(vips[,1]>1))]
write.table(X, file=paste("subset.by.vips", omic, model, "tsv", sep="."), sep="\t", quote = F, row.names = F)

# PCA of the data set with only the selected variables (PERMANOVA done separatedly)
X = t(X)
pca = pca(X, center = F, scale = F)
png(paste0("pca_",omic,".png"))
plotIndiv(pca, group = splsda.res$Y, ind.names = F, pch = c(16,16), title = paste('PCA -', omic, "-", model),
          star = T, legend = T, legend.position = "bottom")
dev.off()

##########
## Significant variables
# t-test and p-value adjustment
pvalues = c()
for (i in 1:ncol(X)) {
  p = t.test(X[rownames(X) %in% meta$pares[which(meta$grupos=="increase")], i], 
             X[rownames(X) %in% meta$pares[which(meta$grupos=="decrease")], i])
  pvalues = c(pvalues, p$p.value)
  names(pvalues)[i] = colnames(X)[i]}
padj = p.adjust(pvalues, method = "BH")
pvalues = data.frame(metabolite = names(pvalues), p.value=pvalues, p.value.adj=padj)
pvalues = pvalues[order(pvalues$p.value.adj),]
write.table(pvalues, file=paste(omic,"pvalues", model, "tsv", sep="."), sep="\t", quote = F, row.names = F)

# Boxplots for significant metabolites
if (omic == "metabolomics" | omic == "16S rRNA analysis") {
  select_for_boxplot = as.character(pvalues$metabolite[which(pvalues$p.value.adj<=0.1 & pvalues$p.value<=0.05)])
  for (i in select_for_boxplot){
    metab = i
    q = round(pvalues[pvalues$metabolite==metab, "p.value.adj"], digits = 4)
    p = round(pvalues[pvalues$metabolite==metab, "p.value"], digits = 4)
    metabolito = X[,colnames(X) == metab]
    if (model == "predictive") metabolito = metabolito + log2(1000000) ## Transform to counts per million
    metabolito=as.data.frame(cbind(valores = metabolito, grupos = meta[names(metabolito), "grupos"]))
    metabolito$valores = as.numeric(as.vector(metabolito$valores))
    metab2 = gsub(pattern = "\\.[0-9]", replacement = "", x=metab)
    metab2 = paste0(gsub(pattern = "\\.", replacement = " ", x=metab2))
    metab2 = stringr::str_to_sentence(metab2)
    if (metab2 == "Phanylalanine") metab2="Phenylalanine"
    if (metab2 == "Ile") metab2 = "Isoleucine"
    if (metab2 == "Val") metab2 = "Valine"
    if (metab2 == "5-aminovaleric acid") metab2 = "5-aminovalerate"

    pl <- ggplot(metabolito, aes(x=factor(grupos), y=valores))
    p2 = pl +
      geom_boxplot(aes(fill=factor(grupos))) +
      theme_bw() +
      labs(fill="Grupos") +
      stat_summary(fun=mean, geom="point", shape=20, size=1, color="red", fill="red") +
      xlab("Groups") +
      ylab("Values per sample") +
      geom_point( aes(x=factor(grupos),y=valores), alpha=0.5 ) +
      ggtitle(paste0(metab2, " (p = ", p, ", q = ",q,")")) +
      theme(plot.title = element_text(size = 25, face = "bold"),
            axis.text.y = element_text(size=20),
            axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            legend.position = "none")+
      scale_fill_manual(values = c("MRE decrease" = "red2",
                                   "MRE increase" = "forestgreen"))
    png(file=paste0("outliers.",metab,".png"))
    print(p2)
    dev.off() 
    
  }
}

# Heatmaps and enrichment for metagenomics and proteomics
if (omic == "metagenomics" | omic = "proteomics"){
  # HEATMAPS
  otus_matrix = t(X[,colnames(X) %in% as.character(pvalues[pvalues$p.value<=0.05 & pvalues$p.value.adj<=0.1,1])])
  colnames(otus_matrix) = make.names(meta[meta$pares %in% colnames(otus_matrix), "grupos"], unique = T)
  otus_matrix = otus_matrix[,order(colnames(otus_matrix))]
  M=otus_matrix
  M=as.matrix(M)
  storage.mode(M) = "numeric"
  data.temp <- M
  for (i in 1:nrow(data.temp)){
    minimo = min(data.temp[i,which(data.temp[i,] != 0)])
    data.temp[i,] = data.temp[i,] + abs(minimo)
    data.temp[i,] = data.temp[i,] / sum(data.temp[i,])
  }
  
  data.temp <- t(scale(t(data.temp)))

  hc1 <- hclust(dist(t(data.temp[,grepl(pattern = "^d.+", x=colnames(data.temp))])), method = "complete")
  hc2 <- hclust(dist(t(data.temp[,grepl(pattern = "^i.+", x=colnames(data.temp))])), method = "complete")
  hc = merge(as.dendrogram(hc1),rev(as.dendrogram(hc2)))
  
  hr <- hclust(dist(data.temp), method = "complete")
  lmat = rbind(c(0,3),c(2,1),c(0,4))
  lwid = c(1.5,4)
  lhei = c(2,4,1)
  colors = ifelse(grepl("decrease.*", colnames(data.temp)), "red2", "forestgreen")
  myBreaks <- seq(-2, 2, length.out=11)
  pdf(paste0(omica, "-", filtrado, log,"-heatmap.version_final",comp,".pdf"))
  heatmap.2(data.temp,breaks = myBreaks,col=colorRampPalette(c("red","yellow","darkgreen"))(10), Colv=as.dendrogram(hc),
            Rowv=as.dendrogram(hr),dendrogram="both",
            trace="none", key=F, keysize = 1.5,lhei = c(2,4), lwid = c(0.5,1),
            ylab = "KO", xlab="Samples", ColSideColors = colors, symkey=FALSE, density.info="density",
            labRow = F, labCol = F, main = paste0("Heatmap - ", omic, " - ", model))
  
  legend(locator(), legend = c("MRE decrease", "MRE increase"), col= c("red2", "forestgreen"), lty= 1.5,
         lwd = 2, cex=1, xpd=TRUE)
  dev.off()
  
  #ENRICHMENT
  kegg_path = read.delim("~/tfm/transcriptomics/processed/ko_andpathways_desglosed.tsv", sep="\t", stringsAsFactors = F,
                         header = T)
  kegg_path = as.data.frame(cbind(desc = kegg_path$path_desc, ko = kegg_path$ko))

  colnames(otus_matrix) = make.names(meta[meta$pares %in% colnames(otus_matrix), "grupos"], unique = T)
  df_pos = otus_matrix[,grepl("increase.*", colnames(otus_matrix))]
  df_neg = otus_matrix[,grepl("decrease.*", colnames(otus_matrix))]
  result = data.frame(KEGG = rep(NA, nrow(otus_matrix)), BOOL = rep(NA,nrow(otus_matrix)))
  
  # Select variables whether they are increased in MRE increase group or not
  for (i in 1:nrow(otus_matrix)){
    result[i,1] = as.character(rownames(otus_matrix)[i])
    result[i,2] = mean(as.numeric(df_pos[i,])) > mean(as.numeric(df_neg[i,]))
  }
  
  info_kegg = result[result$BOOL == "TRUE",1]
  ewp <- enricher(info_kegg, TERM2GENE = kegg_path, pvalueCutoff = 0.1)
  barplot(ewp) + ylab("Number of KOs per KEGG pathway") + ggtitle(paste0("Enrichment analysis - ", omic, " - ", model)) +
    theme(plot.title.position = "plot",
          plot.title = element_text(size=20))
  
  info_kegg = result[result$BOOL == "FALSE",1]
  ewp <- enricher(info_kegg, TERM2GENE = kegg_path, pvalueCutoff = 0.1)
  barplot(ewp) + ylab("Number of KOs per KEGG pathway") +
    ggtitle(ggtitle(paste0("Enrichment analysis - ", omic, " - ", model))) + theme(plot.title.position = "plot",
          plot.title = element_text(size=20))
}



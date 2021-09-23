#!/usr/bin/env Rscript

parameters <- commandArgs(trailingOnly=TRUE)

counts.file <- as.character(parameters[1]) # Omic dataset
meta.file <- as.character(parameters[2]) # Metadata file
variable <- as.character(parameters[3]) # Variable to contrast
file.struc.zeros <- as.character(parameters[4]) # Table with other normalization
                                                  # to do boxplots of structural zeros (CLR+1)
pacman::p_load(dplyr, tidyr, ggplot2)

ancom2<-function(OTUdat, Vardat, main.var="country", pr=NULL, pii=0.5,ref_name=NULL){
  source("pop_detect.R");  source("stepA1.R");  source("stepA2.R");  source("geom_ref.R")
  ########## STEP A1 ################################################################
  A1=struc_zero(OTUdat = OTUdat, Vardat = Vardat, main.var = main.var, p=pr)
  ########## STEP A2 ################################################################
  A2=reset_dat_missing(dat=A1[[1]],pii=pii,ref_name=ref_name,zeros=A1[[2]])
  adjusted_data=A2[[1]]; number_A2_adj_microbes=A2[[2]]; names_A2_adj_microbes=A2[[3]]
  ########## STEP A3 ################################################################
  adjusted_data[adjusted_data==0]=1
  ########## STEP B ################################################################
  ########## filter struc zero################################################################
  struc_zero_elim_ind= which(colSums(A1[[2]])>1);   sid=adjusted_data[,1];  
  normalu=adjusted_data[,3];   popu=adjusted_data[,2];  dat=adjusted_data[,-c(1,2,3)]
  ###Normalize and test individually
  prueba = dat[,-struc_zero_elim_ind]
  prueba = apply(prueba, 2, function(x) log(as.numeric(x)))
  redat = apply(prueba, 2, function(x) x-as.numeric(normalu))
  gauss_data=data.frame(Sample.ID=sid, pop=popu, redat)
  #############subset only the microbes#######################################
  gauss_microbiome=gauss_data %>% select(-c(Sample.ID,pop))
  ########################analysis ##########################################
  d=dim(gauss_microbiome)[2]; pval=matrix(0,d,1); sig_ind=matrix(0,d,1)
  for (j in 1:d){covariate=gauss_data$pop; response=gauss_microbiome[,j]
  regress=data.frame(response,covariate);  regress2=na.omit(regress)
  if (length(unique(regress2$covariate)) > 1){
  model=lm(regress2$response~ factor(regress2$covariate));  anv=anova(model);  pval[j]=anv$`Pr(>F)`[[1]]}}
  padj=p.adjust(pval,"BH"); ind=which(padj<0.1)
  detected_microbes=names(gauss_microbiome)[ind]
  detected_microbes = as.data.frame(cbind(detected_microbes = detected_microbes, pval = pval[ind], padj = padj[ind]))
  structural_zeros= data.frame(Microbe=names(adjusted_data[-c(1,2,3)]),t(1-A1[[2]]))
  gauss_data = t(gauss_data)
  colnames(gauss_data) = gauss_data[1,]
  gauss_data = gauss_data[-c(1,2),]
  gauss_data = cbind(ID = rownames(gauss_data), gauss_data)
  structural_zeros$sum = structural_zeros$yes+structural_zeros$no
  return(list(norm_data=gauss_data, detected_microbes=detected_microbes, structural_zeros=structural_zeros))}
####################################################################################

OTUdat = read.table(counts.file, header = T, sep="\t", stringsAsFactors = F, check.names = F, row.names = 1)
Vardat = read.table(meta.file, header = T, sep="\t", stringsAsFactors = F, check.names = F)
comun = intersect(colnames(OTUdat), Vardat$sample)
OTUdat = OTUdat[,comun]
Vardat = Vardat[Vardat$sample %in% comun,]
names(Vardat)[1] = "Sample.ID"
OTUdat = t(OTUdat)
OTUdat = cbind(Sample.ID = rownames(OTUdat), OTUdat)
final = ancom2(OTUdat, Vardat, main.var=variable)

norm_data = as.data.frame(final$norm_data)
counts.file = unlist(strsplit(x=counts.file, split = "/"))[length(unlist(strsplit(x=counts.file, split="/")))]
write.table(norm_data, file = paste0("ancomii.norm.",counts.file), sep = "\t", quote = F, row.names = F)

detected_microbes = final$detected_microbes
detected_microbes = detected_microbes[order(as.numeric(detected_microbes$pval)),]
detected_microbes$log2fc = sapply(detected_microbes$detected_microbes, function(x){
  log2(mean(as.numeric(norm_data[norm_data$ID==x, colnames(norm_data) %in% Vardat$Sample.ID[which(Vardat[,variable]=="yes")]]))/
        mean(as.numeric(norm_data[norm_data$ID==x, colnames(norm_data) %in% Vardat$Sample.ID[which(Vardat[,variable]=="no")]])))
  })
write.table(detected_microbes, file = paste0("detected_microbes.tsv"), sep = "\t", quote = F, row.names = F)


structural_zeros = final$structural_zeros
write.table(structural_zeros, file = paste0("structural_zeros.tsv"), sep = "\t", quote = F, row.names = F)

select_for_boxplot = as.character(detected_microbes$detected_microbes[which(detected_microbes$pval>0 & detected_microbes$pval <=0.05
                                                                            & detected_microbes$padj>0& detected_microbes$padj <=0.1)])
norm_data = t(norm_data)
norm_data = norm_data[-1,]
pdf(file=paste0("ancomii.plots",".pdf"))
for (i in select_for_boxplot){
  metab = i
  q = round(as.numeric(detected_microbes[detected_microbes$detected_microbes==metab, "padj"]), digits = 4)
  p = round(as.numeric(detected_microbes[detected_microbes$detected_microbes==metab, "pval"]), digits = 4)

  
  metabolito = as.numeric(norm_data[,colnames(norm_data) == metab])
  names(metabolito) = rownames(norm_data)
  metabolito=as.data.frame(cbind(valores = metabolito, grupos = Vardat[Vardat$Sample.ID %in% names(metabolito), variable]))
  metabolito$valores = as.numeric(as.vector(metabolito$valores))
  fc = log2(mean(metabolito$valores[metabolito$grupos=="yes"])/mean(metabolito$valores[metabolito$grupos=="no"]))
  
  if (nchar(metab)>10) metab = unlist(strsplit(metab, split="__"))[length(unlist(strsplit(metab, split="__")))]
  pl <- ggplot(metabolito, aes(x=factor(grupos), y=valores))
  p2 = pl +
    geom_boxplot(aes(fill=factor(grupos))) +
    theme_bw() +
    labs(fill="Grupos") +
    stat_summary(fun=mean, geom="point", shape=20, size=5, color="forestgreen", fill="forestgreen") +
    xlab("Groups") +
    ylab("Values per sample") +
    geom_point( aes(x=factor(grupos),y=valores), alpha=0.5 ) +
    ggtitle(paste0(metab, " (p = ", p, ", q = ",q,")", "\nlog2fc(yes/no) = ",fc)) +
    theme(plot.title = element_text(size = 20, face = "bold"),
          axis.text.y = element_text(size=20),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank())+
    scale_fill_manual(values = c("yes" = "red2",
                                 "no" = "blue"))
  print(p2)
  
}
dev.off() 

###### BOXPLOTS IN CLR NORM OF VARIABLES WITH STRUCTURAL ZEROS ######
table.struc.zeros = read.table(file.struc.zeros, header = T, sep="\t", stringsAsFactors = F, row.names = 1, check.names = F)
comun = intersect(colnames(table.struc.zeros), Vardat$Sample.ID)
table.struc.zeros = table.struc.zeros[,comun]
table.struc.zeros = as.data.frame(t(table.struc.zeros))
colnames(table.struc.zeros) = gsub(x=colnames(table.struc.zeros), pattern = "-", replacement = ".")
colnames(table.struc.zeros) = gsub(x = colnames((table.struc.zeros)), pattern = ";", replacement = ".")
select_for_boxplot = intersect(detected_microbes$detected_microbes, colnames(table.struc.zeros))
pdf(file=paste0("clr.plots",".pdf"))
for (i in select_for_boxplot){
  metab = i
  metabolito = as.numeric(table.struc.zeros[,colnames(table.struc.zeros) == metab])
  names(metabolito) = rownames(table.struc.zeros)
  metabolito=as.data.frame(cbind(valores = metabolito, grupos = Vardat[Vardat$Sample.ID %in% names(metabolito), variable]))
  metabolito$valores = as.numeric(as.vector(metabolito$valores))
  if (nchar(metab)>10) metab = unlist(strsplit(metab, split="__"))[length(unlist(strsplit(metab, split="__")))]
  pl <- ggplot(metabolito, aes(x=factor(grupos), y=valores))
  p2 = pl +
    geom_boxplot(aes(fill=factor(grupos))) +
    theme_bw() +
    labs(fill="Grupos") +
    stat_summary(fun=mean, geom="point", shape=20, size=5, color="forestgreen", fill="forestgreen") +
    xlab("Groups") +
    ylab("Values per sample") +
    geom_point( aes(x=factor(grupos),y=valores), alpha=0.5 ) +
    ggtitle(paste0(metab)) +
    theme(plot.title = element_text(size = 20, face = "bold"),
          axis.text.y = element_text(size=20),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank())+
    scale_fill_manual(values = c("yes" = "red2",
                                 "no" = "blue"))
  print(p2)
  
}
dev.off() 

#' @title Eset50467
#' 
#' @description
#' Effect of SRSF2 over-expression in a human lung carcinoma cell line
#' 
#' Homo sapiens.
#' 
#' Expression Set from experiment GSE50467: Effect of SRSF2 over-expression in a human lung carcinoma cell line.
#'   SRSF2 is an important gene for the splicing of a great amount of genes.
#' Additionally, as it can bee seen in this study, the overexpression of this gene change the expression of a lot
#'   of genes.
#' The overall design includes 6 samples with SRSF2 over-expression and 6 samples as control of H358 cell line.
#' Expression data has been background-corrected and normalized.
#' In addition to expression data, phenoData, experimentData and featureData can also be found.
#' @source
#' \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50467}
#' @examples
#' data(Eset50467,package="enriquepresa")
#' @docType data
#' @keywords datasets
#' @format Large ExpressionSet
#' 
#' @name Eset50467
NULL


#' @title genefilter_df
#' 
#' @description
#' Dataframe generated with the GSE50467 expression set.
#' This dataframe includes differential expression information for each probe.
#' The analysis was done with a t-test with equal variances assumed using the function rowttests
#'   of genefilter package.
#' Dataframe includes the gene symbol, the t statistic, the LogFC, the p-value, the adjusted
#'   p-value with both Benjamini-Hochberg and Bonferroni method, and the URL's to ENTREZ and ENSEMBL.
#' @source
#' \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50467}
#' @examples
#' data(genefilter_df,package="enriquepresa")
#' @docType data
#' @keywords datasets
#' @format data.frame
#' 
#' @name genefilter_df
NULL


#' @title limma_df
#' 
#' @description
#' Dataframe generated with the GSE50467 expression set.
#' This dataframe includes differential expression information for each probe.
#' The analysis was done with a t-test without equal variances assumed using limma package.
#' Dataframe includes the gene symbol, the t statistic, the LogFC, the p-value the adjusted p-value
#'   with both Benjamini-Hochberg and Bonferroni method, and the URL's to ENTREZ and ENSEMBL.
#' @source
#' \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50467}
#' @examples
#' data(limma_df,package="enriquepresa")
#' @docType data
#' @keywords datasets
#' @format data.frame
#' 
#' @name limma_df
NULL


#' @title PRJNA517180
#' 
#' @description
#' A Myt1 family transcription factor defines neuronal fate by repressing non-neuronal genes.
#' 
#' Caenorhabditis elegans.
#' 
#' RNA-seq data from PRJNA51780. In this experiment, the importance of ztf-11 for neurogenesis was evaluated in
#'   order to infer the role of Myt1, the Homo sapiens homologous gene.
#' Myt1 and ztf-11 are zinc-finger transcription factors that are thought to be important in the reprogrammation
#'   from fibroblast to neuron in vitro.
#' There are two groups: the control group consists of individuals with a control vector (L4440), which does not
#'   affect physiology.
#' On the other hand, the case group consists of knockdown mutants of ztf-11 gene via RNAi.
#' @source
#' \url{https://www.ncbi.nlm.nih.gov/bioproject/PRJNA517180}
#' @examples
#' data(PRJNA51780,package="enriquepresa")
#' @docType data
#' @keywords datasets
#' @format Large RangedSummarizedExperiment
#' 
#' @name PRJNA517180
NULL

#' @title edgeR_df
#' 
#' @description
#' Data frame with the result of doing a differential expression analysis with edgeR package and the
#'    SummarizedExperiment PRJNA517180. This data frame is included as data in order to be able to do the
#'    overrepresentation analysis. 
#' In this data frame, the WORMBASE identifiers as rownames, gene symbols, logFC, logCPM, P-value,
#'    Benjamini-Hochberg adjusted P-value, and ENTREZ and ENSEMBL URLs are included.
#' @source
#' \url{https://www.ncbi.nlm.nih.gov/bioproject/PRJNA517180}
#' @examples
#' data(edgeR_df,package="enriquepresa")
#' @docType data
#' @keywords datasets
#' @format data.frame
#' 
#' @name edgeR_df
NULL

#' @title ENTREZ URL generation
#' URL for an ENTREZ ID
#' @description  
#' It makes the url from ENTREZ giving the ENTREZ ID.
#' @param entrezid ENTREZ ID for a gene
#' @return 
#' A character referring to the ENTREZ URL
#' @export
#' @family URL generation
entrezurl = function(entrezid){
  ifelse(entrezid == "NA", NA, paste0("<a href='http://www.ncbi.nlm.nih.gov/gene/?term=", 
                               entrezid, "'>", entrezid, "</a>"))
}

#' @title ENSEMBL URL generation
#' URL for an ENSEMBL ID
#' @description  
#' It makes the url from ENSEMBL giving the ENSEMBL ID.
#' @param ensemblid ENSEMBL ID for a gene
#' @return 
#' A character referring to the ENSEMBL URL
#' @export
#' @family URL generation
ensemblurl = function(ensemblid){
  ifelse(ensemblid == "NA", NA, paste0("<a href='http://www.ncbi.nlm.nih.gov/gene/?term=", 
                               ensemblid, "'>", ensemblid, "</a>"))
}

#' @title generation of a data frame with URLs from GO
#' URLs for a set of GO terms
#' @description  
#' This function receives a data frame that is the result of some GO enrichment analysis and makes a
#'   data frame with URLs from GO for every GO term. The GO terms must be as rownames.
#' @param df data frame that is the result of some GO enrichment analysis and with GO terms as rownames
#' @return 
#' A data frame with URLs from GO for every GO term
#' @export
#' @family URL generation
dfGO = function(df){
  lista = sapply(rownames(df), function(x){expr = regexpr("GO:[0-9]+",x)
  return(regmatches(x, expr))})
  
  GO.URLS = sapply(lista, function(x) return(paste0("<a href='http://amigo.geneontology.org/amigo/term/", x, "'>", x, "</a>")))
  return(cbind(df, GO.URLS))
}
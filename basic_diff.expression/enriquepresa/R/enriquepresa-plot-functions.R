#' @title Volcano plot generation
#' @description  
#' It generates a Volcano plot for a set of genes. A volcano plot is a kind of representation where the genes
#'   that are differentially expressed over a p-value (here represented as the opposite decimal logarithm of the
#'   p-value) are remarked. 
#' @param adj_p information of adjusted p-values of a dataframe
#' @param logFC logFC for each gene of the same dataframe
#' @param symbols gene names in order to mark those that are differentially expressed 
#' @param min_p adjusted p-value that is the limit of significance. By default is 0.05
#' @param diff limit of logFC to consider a gene interesting. By default is 1 (a gene must be two times
#'   overexpressed or underexpressed)
#' @return 
#' Volcano plot 
#' @export
#' @family Volcano plot generation
plot_volcano = function(adj_p, logFC, symbols, min_p=0.05, diff=1){
  logs = ifelse(adj_p == 0, 0, -log10(adj_p))
  plot(logFC,logs,
       xlim=range(logFC, na.rm=TRUE),
       ylim=range(logs, na.rm=TRUE),
       xlab="LogFC",
       ylab="-log10(p-value)")
  abline(h=-log10(min_p),col="green")
  abline(v=diff,col="blue")
  abline(v=-diff,col="blue")
  seleccion = which(abs(logFC) > diff & logs > -log10(min_p))
  points(logFC[seleccion],logs[seleccion],pch=18,col="red")
  text(logFC[seleccion],logs[seleccion], symbols[seleccion], cex = 0.5, pos = 1)
  title("Volcano plot")
}

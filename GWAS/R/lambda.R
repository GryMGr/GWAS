#' Calculate genomic inflation factor for GWAS/EWAS resutls.
#'
#' @param x A GWAS or EWAS result data matrix, with the p-values as the third
#' column. If the p-values are not in 3rd col, than col has to be set to what-
#' ever column the p-values are located in. 
#' @param col Column indicating where the p-values are located. Default is set
#' to 3.
#' @return Returns the genomic inflation factor of an set of genomic tests. 
#' Caution should be made when interpreting the lambdas of epigenetic results
#' since they not always reflects the biology or epidemiology in question for 
#' epigenetic analysis. 
#' @export
#' @example
#' lambda(x)
lambda <- function(x, col = 3){
	return(qchisq(median(x[,col]), df = 1, lower.tail = FALSE)/qchisq(0.5, df = 1))
}

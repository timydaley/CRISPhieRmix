#' Simulated data from Rosenbluh et al (2017) Nature Communications
#'
#' A dataset containing DESeq2 computed log2 fold changes for 34 essential genes
#' and 102 simulated negative genes, as well as negative control guides.
#' Simulated negative genes have geneIds starting with 'sim'
#'
#' @name Rosenbluh2017CRISPRiSim
#' @docType data
#' @keywords data
#' @format a list with log2fc = log2 fold changes of gene targetting guides, 
#'   geneIds = gene ids of gudies,
#'   negCtrl = log2 fold changes of negative control guides, 
#'   and counts = the count matrix 
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5457492/}
#'
#' @examples
#' data(Rosenbluh2017CRISPRiSim)
#' 
"Rosenbluh2017CRISPRiSim"

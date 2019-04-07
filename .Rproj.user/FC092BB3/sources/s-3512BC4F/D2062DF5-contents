#'The process of geting the differential expression genes
#'
#' This function uses based method such as limma to get the differential expression gene
#'
#' the previous SPIA method and integrate the change of
#' of genes Pearson coefficient(PCC) from two groups. We proposed a set of three
#' pathway analysis methods based on the change of PCC. We applied these approaches
#' to colorectal cancer, lung cancer and Alzheimer's disease datasets and so on.
#'
#' @param dateset   The dataset of cancer or disease from NCBI, such as GSE1145
#' @export
#' @return There are five variable values, exprs stand for gene expression profile,
#' normal stand for the number of normal samples, tumor stand for the number of tumor
#' samples, DE stand for the fold change of differential expression genes, ALL stand for
#' all genes of human
#' @examples
#' #import EnrichmentBrowser, KEGGandMetacoreDzPathwaysGEO, KEGGdzPathwaysGEO and SPIA package
#' library(EnrichmentBrowser)
#' library(KEGGandMetacoreDzPathwaysGEO)
#' library(KEGGdzPathwaysGEO)
#' library(SPIA)
#' data("GSE1145")
#' result <- process(GSE1145)


process <- function(dataset){

  #change probe to gene symbol
  all.eset <- probe.2.gene.eset(dataset)
  # Normalization of gene expression profile

  all.eset <- normalize(all.eset, norm.method="quantile")
  after.norm <- exprs(all.eset)

  exprs_norm <- data.frame(after.norm)
  #table(pData(all.eset)$Group)
  pData(all.eset)$GROUP <- ifelse(pData(all.eset)$Group == "d", 1, 0)
  normal <- length(which(pData(all.eset)$GROUP == '0'))
  tumor <- length(which(pData(all.eset)$GROUP == '1'))

  # get differential expression genes
  all.eset <- de.ana(all.eset)
  head(fData(all.eset), n=4)
  all_de <- fData(all.eset)

  tg <- all_de[all_de$ADJ.PVAL < 0.1,]
  # get fold change pf differential expression genes
  DE = tg$FC
  names(DE)<-as.vector(rownames(tg))
  # get all gene names
  ALL = rownames(all_de)

  res = list(exprs = exprs_norm, normal = normal, tumor = tumor,
             DE = DE, ALL = ALL)
  return(res)

}



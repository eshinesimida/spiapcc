#' Calculate spiapcc pathway scores from gene expression base on the change of PCC
#'
#' This function uses the previous SPIA method and integrate the change of
#' of genes Pearson coefficient(PCC) from two groups. We proposed a set of three
#' pathway analysis methods based on the change of PCC. We applied these approaches
#' to colorectal cancer, lung cancer and Alzheimer's disease datasets and so on.
#'
#'
#' We used a compendium of 22 GEO datasets obtained from the KEGGdzPathwaysGEO
#' and the KEGGandMetacoreDzPathwaysGEO benchmark sets(Tarca et al., 2012; Tarca et al., 2013)
#' in this study. These datasets have been specifically chosen as they study a
#' certain human disease in a corresponding KEGG pathway (e.g. Alzheimer's disease).
#' These pathways are regarded as the target pathways in the following. We investigate
#' first how well the individual set- and network-based methods detect the target
#'
#' The expression data sets involved 11 conditions and 17 tissues.
#'
#'
#' @param de   The number of differential genes
#' @param all  All genes in human
#' @param normal the number of normal samples
#' @param tumor the number of tumor samples
#' @param flag flag = 1,0,-1 , if flag = 1 from normal to tumor, flag = -1 from tumor to normal, flag = 0 stand for absolute value between tow groups
#' @export
#' @examples
#' #import EnrichmentBrowser, KEGGandMetacoreDzPathwaysGEO, KEGGdzPathwaysGEO and SPIA package
#' library(EnrichmentBrowser)
#' library(KEGGandMetacoreDzPathwaysGEO)
#' library(SPIA)
#' data("GSE8671")
#' # Get expression profile of GSE8671
#' exprs_all <- exprs(GSE8671)
#' # Add the gene symbol
#' all.eset <- probe.2.gene.eset(GSE8671)
#' head(featureNames(all.eset))
#' before.norm <- exprs(all.eset)
#' # Gene normalization
#' all.eset <- normalize(all.eset, norm.method="quantile")
#' after.norm <- exprs(all.eset)
#' exprs_all1 <- data.frame(after.norm)
#' table(pData(all.eset)$Group)
#' pData(all.eset)$GROUP <- ifelse(pData(all.eset)$Group == "d", 1, 0)
#' normal <- length(which(pData(all.eset)$GROUP == '0'))
#' tumor <- length(which(pData(all.eset)$GROUP == '1'))
#' # Get the differential expression genes in limmar package
#' all.eset <- de.ana(all.eset)
#' head(fData(all.eset), n=4)
#' all_de <- fData(all.eset)
#' tg <- all_de[all_de$ADJ.PVAL < 0.1,]
#' DE_colorectal = tg$FC
#' names(DE_colorectal)<-as.vector(rownames(tg))
#' ALL_colorectal = rownames(all_de)
#' #The result of spia method
#' res_spia = spia(de = DE_colorectal, all=ALL_colorectal, organism="hsa",nB=2000,plots=FALSE,beta=NULL,combine="fisher",verbose=TRUE)
#' gse_madat2 <- exprs_all1
#' # The results of spia_nt method
#' res_nt = spiap(de=DE_colorectal, all=ALL_colorectal, gse_madat2 = gse_madat2,normal = normal, tumor = tumor,norganism="hsa",nB=2000,plots=FALSE,
#'                   beta=NULL,combine="fisher",verbose=T, flag = 1)
#' #The results of spia_tn method
#' res_tn = spiap(de=DE_colorectal, all=ALL_colorectal, gse_madat2 = gse_madat2, normal = normal, tumor = tumor,organism="hsa",nB=2000,plots=FALSE,
#'                   beta=NULL,combine="fisher",verbose=T, flag = -1)
#' #The results of spia_abs method
#' res_abs = spiap(de = DE_colorectal, all=ALL_colorectal , gse_madat2 = gse_madat2,normal = normal, tumor = tumor,organism="hsa",nB=2000,plots=FALSE,
#'                   beta=NULL,combine="fisher",verbose=T, flag = 0)

#spiap <- function(de = NULL, all = NULL,gse_madat2 = gse_madat2, normal = NULL, tumor = NULL, organism = "hsa", data.dir = NULL,
#                  pathids = NULL, nB = 2000, plots = FALSE, verbose = TRUE,
#                  beta = NULL, combine = "fisher"){
 # UseMethod("spiap")
#}

spiapcc <- function (de = NULL, all = NULL,gse_madat2 = gse_madat2, normal = NULL, tumor = NULL, organism = "hsa", data.dir = NULL,
                       pathids = NULL, nB = 2000, plots = FALSE, verbose = TRUE,
                       beta = NULL, combine = "fisher", flag = -1)
{
  if (is.null(de) | is.null(all)) {
    stop("de and all arguments can not be NULL!")
  }
  rel <- c("activation", "compound", "binding/association",
           "expression", "inhibition", "activation_phosphorylation",
           "phosphorylation", "inhibition_phosphorylation", "inhibition_dephosphorylation",
           "dissociation", "dephosphorylation", "activation_dephosphorylation",
           "state change", "activation_indirect effect", "inhibition_ubiquination",
           "ubiquination", "expression_indirect effect", "inhibition_indirect effect",
           "repression", "dissociation_phosphorylation", "indirect effect_phosphorylation",
           "activation_binding/association", "indirect effect",
           "activation_compound", "activation_ubiquination")
  if (is.null(beta)) {
    beta = c(1, 0, 0, 1, -1, 1, 0, -1, -1, 0, 0, 1, 0, 1,
             -1, 0, 1, -1, -1, 0, 0, 1, 0, 1, 1)
    names(beta) <- rel
  }
  else {
    if (!all(names(beta) %in% rel) | length(names(beta)) !=
        length(rel)) {
      stop(paste("beta must be a numeric vector of length",
                 length(rel), "with the following names:", "\n",
                 paste(rel, collapse = ",")))
    }
  }
  .myDataEnv <- new.env(parent = emptyenv())
  datload <- paste(organism, "SPIA", sep = "")
  if (is.null(data.dir)) {

    load(file = paste(system.file("extdata", package = "SPIA"),
                      paste("/", organism, "SPIA", sep = ""), ".RData",
                      sep = ""), envir = .myDataEnv)
  }

  datpT = .myDataEnv[["path.info"]]

  datp <- list()
  path.names <- NULL
  hasR <- NULL

  for (jj in 1:length(datpT)) {
    sizem <- dim(datpT[[jj]]$activation)[1]
    s <- 0
    con <- 0
    for (bb in 1:length(rel)) {
      con = con + datpT[[jj]][[rel[bb]]] * abs(sign(beta[rel[bb]]))
      s = s + datpT[[jj]][[rel[bb]]] * beta[rel[bb]]
    }
    ####################
    ####Wij produce ####
    ####################
    id <- rownames(s)
    gene_exp1 <- gse_madat2[match(id, rownames(gse_madat2)),1:normal]
    gene_exp2 <- gse_madat2[match(id, rownames(gse_madat2)),(normal+1):(normal+tumor)]
    gene_exp1_mat <- cor(t(gene_exp1))
    diag(gene_exp1_mat) <- 0

    gene_exp2_mat <- cor(t(gene_exp2))
    diag(gene_exp2_mat) <- 0
    # difference from tumor to normal
    if(flag == -1){
      diff_exp <- -(gene_exp1_mat - gene_exp2_mat)
    }
    # difference form normal to tumor
    else if(flag == 1){
      diff_exp <- -(gene_exp1_mat - gene_exp2_mat)
    }
    # absolute value of difference between tow groups
    else if(flag == 0){
      diff_exp <- abs(gene_exp1_mat - gene_exp2_mat)
    }


    id_1 <- which(s != 0)
    for(i in 1:length(id_1)){
      {
        if(id_1[i] %% dim(s)[1] == 0){
          id2 <- id_1[i]/dim(s)[1]
          id1 <- dim(s)[1]
        }else{
          id2 <- floor(id_1[i]/dim(s)[1] + 1)
          id1 <- id_1[i] %% dim(s)[1]
        }
      }
      if(!is.na(diff_exp[id1,id2])){
        s[id1,id2] <- diff_exp[id1,id2] * s[id1,id2]
      }else{
        s[id1,id2] <- s[id1,id2]
      }
    }

    ########################################################
    #############################
    z = matrix(rep(apply(con, 2, sum), dim(con)[1]), dim(con)[1],
               dim(con)[1], byrow = TRUE)
    z[z == 0] <- 1
    datp[[jj]] <- s/z
    path.names <- c(path.names, datpT[[jj]]$title)
    hasR <- c(hasR, datpT[[jj]]$NumberOfReactions >= 1)
  }
  names(datp) <- names(datpT)
  names(path.names) <- names(datpT)
  tor <- lapply(datp, function(d) {
    sum(abs(d))
  }) == 0 | hasR | is.na(path.names)
  datp <- datp[!tor]
  path.names <- path.names[!tor]
  IDsNotP <- names(de)[!names(de) %in% all]

  ph <- pb <- pcomb <- nGP <- pSize <- smPFS <- tA <- tAraw <- KEGGLINK <- NULL
  set.seed(1)
  if (plots) {
    pdf("SPIAPerturbationPlots.pdf")
  }
  for (i in 1:length(names(datp))) {
    path <- names(datp)[i]
    M <- datp[[path]]
    diag(M) <- diag(M) - 1
    X <- de[rownames(M)]
    noMy <- sum(!is.na(X))
    nGP[i] <- noMy
    okg <- intersect(rownames(M), all)
    ok <- rownames(M) %in% all
    pSize[i] <- length(okg)
    if ((noMy) > 0 & (abs(det(M)) > 1e-07)) {
      gnns <- paste(names(X)[!is.na(X)], collapse = "+")
      KEGGLINK[i] <- paste("http://www.genome.jp/dbget-bin/show_pathway?",
                           organism, names(datp)[i], "+", gnns, sep = "")
      X[is.na(X)] <- 0
      pfs <- solve(M, -X)
      smPFS[i] <- sum(pfs - X)
      tAraw[i] <- smPFS[i]
      if (plots) {
        par(mfrow = c(1, 2))
        plot(X, pfs - X, main = paste("pathway ID=",
                                      names(datp)[i], sep = ""), xlab = "Log2 FC",
             ylab = "Perturbation accumulation (Acc)", cex.main = 0.8,
             cex.lab = 1.2)
        abline(h = 0, lwd = 2, col = "darkgrey")
        abline(v = 0, lwd = 2, col = "darkgrey")
        points(X[abs(X) > 0 & X == pfs], pfs[abs(X) >
                                               0 & X == pfs] - X[abs(X) > 0 & X == pfs], col = "blue",
               pch = 19, cex = 1.4)
        points(X[abs(X) > 0 & X != pfs], pfs[abs(X) >
                                               0 & X != pfs] - X[abs(X) > 0 & X != pfs], col = "red",
               pch = 19, cex = 1.4)
        points(X[abs(X) == 0 & X == pfs], pfs[abs(X) ==
                                                0 & X == pfs] - X[abs(X) == 0 & X == pfs],
               col = "black", pch = 19, cex = 1.4)
        points(X[abs(X) == 0 & X != pfs], pfs[abs(X) ==
                                                0 & X != pfs] - X[abs(X) == 0 & X != pfs],
               col = "green", pch = 19, cex = 1.4)
      }
      ph[i] <- phyper(q = noMy - 1, m = pSize[i], n = length(all) -
                        pSize[i], k = length(de), lower.tail = FALSE)
      pfstmp <- NULL
      for (k in 1:nB) {
        x <- rep(0, length(X))
        names(x) <- rownames(M)
        x[ok][sample(1:sum(ok), noMy)] <- as.vector(sample(de,
                                                           noMy))
        tt <- solve(M, -x)
        pfstmp <- c(pfstmp, sum(tt - x))
      }
      mnn <- median(pfstmp)
      pfstmp <- pfstmp - mnn
      ob <- smPFS[i] - mnn
      tA[i] <- ob
      if (ob > 0) {
        pb[i] <- sum(pfstmp >= ob)/length(pfstmp) * 2
        if (pb[i] <= 0) {
          pb[i] <- 1/nB/100
        }
        if (pb[i] > 1) {
          pb[i] <- 1
        }
      }
      if (ob < 0) {
        pb[i] <- sum(pfstmp <= ob)/length(pfstmp) * 2
        if (pb[i] <= 0) {
          pb[i] <- 1/nB/100
        }
        if (pb[i] > 1) {
          pb[i] <- 1
        }
      }
      if (ob == 0) {
        if (all(pfstmp == 0)) {
          pb[i] <- NA
        }
        else {
          pb[i] <- 1
        }
      }
      if (plots) {
        bwidth = sd(pfstmp)/4
        if (bwidth > 0) {
          plot(density(pfstmp, bw = bwidth), cex.lab = 1.2,
               col = "black", lwd = 2, main = paste("pathway ID=",
                                                    names(datp)[i], "  P PERT=", round(pb[i],
                                                                                       5), sep = ""), xlim = c(min(c(tA[i] -
                                                                                                                       0.5, pfstmp)), max(c(tA[i] + 0.5, pfstmp))),
               cex.main = 0.8, xlab = "Total Perturbation Accumulation (TA)")
        }
        else {
          pfsTab = table(pfstmp)
          plot(as.numeric(names(pfsTab)), as.numeric(pfsTab),
               cex.lab = 1.2, col = "black", main = paste("pathway ID=",
                                                          names(datp)[i], "  P PERT=", round(pb[i],
                                                                                             5), sep = ""), xlim = c(min(c(tA[i] -
                                                                                                                             0.5, pfstmp)), max(c(tA[i] + 0.5, pfstmp))),
               cex.main = 0.8, xlab = "Total Perturbation Accumulation (TA)",
               ylab = "frequency")
        }
        abline(v = 0, col = "grey", lwd = 2)
        abline(v = tA[i], col = "red", lwd = 3)
      }
      pcomb[i] <- combfunc(pb[i], ph[i], combine)
    }
    else {
      pb[i] <- ph[i] <- smPFS[i] <- pcomb[i] <- tAraw[i] <- tA[i] <- KEGGLINK[i] <- NA
    }
    if (verbose) {
      cat("\n")
      cat(paste("Done pathway ", i, " : ", substr(path.names[names(datp)[i]],
                                                  1, 30), "..", sep = ""))
    }
  }
  if (plots) {
    par(mfrow = c(1, 1))
    dev.off()
  }
  pcombFDR = p.adjust(pcomb, "fdr")
  phFdr = p.adjust(ph, "fdr")
  pcombfwer = p.adjust(pcomb, "bonferroni")
  Name = path.names[names(datp)]
  Status = ifelse(tA > 0, "Activated", "Inhibited")
  res <- data.frame(Name, ID = names(datp), pSize, NDE = nGP,
                    pNDE = ph, tA, pPERT = pb, pG = pcomb, pGFdr = pcombFDR,
                    pGFWER = pcombfwer, Status, KEGGLINK, stringsAsFactors = FALSE)
  res <- res[!is.na(res$pNDE), ]
  res <- res[order(res$pG), ]
  rownames(res) <- NULL
  res
}







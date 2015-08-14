

#' Function to create a PCA analysis
#'
#' This function performs a PCA analysis and returns a list with the PCAresults, PCAscores, PCA plot and scree plot
#' @param fluidSCproc fluidSCproc S3 object
#' @param based_on_values values to use, defaults to "log2Ex"
#' @param PCAscaled boolean to scale data prior to PCA or not, defaults to true
#' @param nr_gene_contribution number of top contributing genes in both directions to show on the PCA plot
#' @param nrClust number of clusters you want
#' @param clustMethod choose method of clustering
#' @param colorClust colors for the different cluster groups, number of colors and clusters must be the same
#' @param firstPC first principal component you want to visualize
#' @param secondPC second principal component you want to visualize
#' @param distMethod choose method to create distance matrix
#' @param corrMethod choose correlation algorithm if clustMethod = "Correlation"
#' @param hclustMethod choose clustering method of distance matrix, defaults to average
#' @param NAvalues choose to remove or replace assays with NA values
#' @param centerpoint boolean to show center of the different groups, defaults to TRUE
#' @param segments boolean to show segments on plot, defaults to TRUE
#' @param corcirc boolean to show correlation circle, defaults to FALSE
#' @param corGenes boolean to show genes that corrleate with sample position, defaults to TRUE
#' @param labeling boolean to show labels of samples, defaults to FALSE
#' @param corr_genes_selection selection of genes to see on the plot, defaults to NULL = all genes if corGenes = T
#' @return returns a list with PCAdata, PCAscores, pcaplot and screeplot
#' @export
#' @details NA
#' @examples
#' PCA_SC()

PCA_SC <- function(fluidSCproc,  based_on_values = "log2Ex", PCAscaled = T, nr_gene_contributions = 5,
                   nrClust = 1, clustMethod = c("Kmeans","Hierarchical", "Correlation"), colorClust = "red", firstPC = "PC1", secondPC = "PC2",
                   distMethod = "euclidean", corrMethod = "pearson", hclustMethod = "average", NAvalues = c("remove_assays","replace_with_0"),
                   centerpoint = T, segments = T, corcirc = F, corGenes = T, labeling = F,  corr_genes_selection = NULL) {
  
  
  # load libraries
  library(gridExtra)
  library(grid)
  library(ggplot2)
  
  normFluidCt <- fluidSCproc$data
  normMatrix <- dcast(normFluidCt, formula = Samples ~ Assays, value.var = based_on_values); rownames(normMatrix) <- normMatrix$Samples; normMatrix <- normMatrix[, -1]
  
  ## for merged fluidSCproc objects, what to do with NA values?
  NAvalues <- match.arg(NAvalues)
  
  # remove genes with NA values
  if(NAvalues == "remove_assays") {
    NAgenes <- colnames(normMatrix)[apply(normMatrix, 2, FUN = function(x) any(is.na(x)))]
    normMatrix <- normMatrix[,-(which(colnames(normMatrix) %in% NAgenes))]
  } 
  # replace NA values with 0
  else if(NAvalues == "replace_with_0") {
    normMatrix[is.na(normMatrix)] <- 0
  } else stop("Choose one of the 2 options to take care of NA values")
  
  
  # do Principal Component Analysis
  PCAdata <- prcomp(normMatrix, scale. = PCAscaled)
  
  # geneloadings represent contribution of variables (gene expressions) to individual principal components
  geneloadings <- PCAdata$rotation
  # find the top genes that contribute the most to the PCs
  top_contributors_PC <- function(PC_loadings, nr_genes = nr_gene_contributions) {
    
    result_vec <- NULL
    
    for(PC in colnames(PC_loadings)){
      
      # relative percentage of gene contribution
      contribution_geneloadings <- abs(PC_loadings[,PC,drop = F])/sum(abs(PC_loadings[,PC, drop = F])) * 100
      tempdfr <- contribution_geneloadings[order(PC_loadings[,PC]),,drop=F]
      
      # genes with negative and positive contribution
      neggenes <- paste0(rownames(tempdfr[1:nr_genes,,drop=F])," ", format(tempdfr[1:nr_genes,], digits = 2),"%",collapse = "  ")
      posgenes <- paste0(rownames(tempdfr[nrow(tempdfr): (nrow(tempdfr)-(nr_genes-1)), ,drop=F])," ", format(tempdfr[nrow(tempdfr) : (nrow(tempdfr) - (nr_genes-1)),], digits = 2),"%", collapse = "  ")
      
      res <- paste0(neggenes, "   VS   ", posgenes)
      
      result_vec <- c(result_vec, res)
    }
    
    return(result_vec)
    
  }
  
  contr <- top_contributors_PC(geneloadings)
  
  
  ### SCREE PLOT ###
  # amount of variance explained by the different prinicipal components 
  
  varExpl <- PCAdata$sdev^2/sum(PCAdata$sdev^2)
  
  varExplDfr <- as.data.frame(cbind(PComp = colnames(PCAdata$x), varExpl = varExpl))
  varExplDfr$PComp <- factor(varExplDfr$PComp, levels = varExplDfr$PComp)
  varExplDfr$varExpl <- as.numeric(as.character(varExplDfr$varExpl))
  
  varExplDfr$geneContr <- contr
  
  # screeplot
  pl <- ggplot(varExplDfr[1:10,], aes(x = PComp, y = varExpl))
  pl <- pl + geom_bar(stat = "identity")
  pl <- pl + theme_bw() + coord_flip()
  pl <- pl + labs(list(x = "Principal components",y = "% of variance explained"))
  
  # contribution of genes
  textpl <- ggplot(varExplDfr[1:10,], aes(x = 1, y = PComp, label = geneContr))
  textpl <- textpl + geom_text(size = 3)
  textpl <- textpl + theme_classic()
  textpl <- textpl + labs(list(x="", y = "Principal components"))
  textpl <- textpl + theme(axis.text.x = element_blank())
  
  
  screeplot <- arrangeGrob(pl, textpl, ncol = 2)
  
  
  ### CLUSTER PCA transformed data ###
  if (nrClust != length(colorClust)) 
    stop("you must provide the same number of colors as number of clusters")
  
  scores <- as.data.frame(PCAdata$x)
  
  clustMethod <- match.arg(clustMethod)
  if (clustMethod == "Kmeans") {
    K <- kmeans(scores[, c(firstPC, secondPC)], centers = nrClust, 
                nstart = 25, iter.max = 1000, algorithm = "Hartigan-Wong")
    I <- K$cluster
  }
  
  else if (clustMethod == "Hierarchical") {
    H <- hclust(dist(scores[, c(firstPC, secondPC)], method = distMethod), 
                method = hclustMethod)
    I <- cutree(H, k = nrClust)
  }
  
  else if (clustMethod == "Correlation") {
    distfun <- function(x) as.dist(1 - cor(t(x), method = corrMethod))
    mycordist <- distfun(scores[, c(firstPC, secondPC)])
    hclustfun <- function(x) hclust(x, method = hclustMethod)
    Hcor <- hclustfun(mycordist)
    I <- cutree(Hcor, k = nrClust)
  }
  
  scores$clust <- I
  
  
  ## Calculate center of groups and segment coordinates ##
  xvecl <- list()
  yvecl <- list()
  for (group in unique(scores$clust)) {
    x <- mean(scores[, firstPC][which(scores$clust == group)])
    y <- mean(scores[, secondPC][which(scores$clust == group)])
    xvecl[paste0("group", group)] <- x
    yvecl[paste0("group", group)] <- y
  }
  scores$xend <- unlist(lapply(scores$clust, function(x) xvecl[paste0("group", x)]))
  scores$yend <- unlist(lapply(scores$clust, function(x) yvecl[paste0("group", x)]))
  
  
  ## variance explained per PC and top contribution of genes ##
  totalsdev <- PCAdata$sdev^2/sum(PCAdata$sdev^2)
  names(totalsdev) <- colnames(PCAdata$rotation)
  
  firstPCgenes <- varExplDfr[varExplDfr$PComp == firstPC,]$geneContr
  secondPCgenes <- varExplDfr[varExplDfr$PComp == secondPC,]$geneContr
  
  firstPCsdev <- format((as.vector(totalsdev[firstPC]) * 100),digits = 2)
  secondPCsdev <- format((as.vector(totalsdev[secondPC]) * 100), digits = 2)
  
  Xtitle <- paste0(firstPC, ": ", firstPCsdev, "% variance")
  Xsubtitle <- paste0("Top genes: ", firstPCgenes)
  Ytitle <- paste0(secondPC, ": ", secondPCsdev, "% variance")
  Ysubtitle <- paste0("Top genes: ", secondPCgenes)
  
  ## circle with correlations ##
  opt_r <- min(max(scores[, firstPC]), max(scores[, secondPC]))
  circle <- function(center = c(0, 0), npoints = 100) {
    r = opt_r
    tt = seq(0, 2 * pi, length = npoints)
    xx = center[1] + r * cos(tt)
    yy = center[1] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
  }
  corcir = circle(c(0, 0), npoints = 100)
  correlations <- as.data.frame(cor(normMatrix, PCAdata$x))
  correlations[, firstPC] <- correlations[, firstPC] * opt_r
  correlations[, secondPC] <- correlations[, secondPC] * opt_r
  
  # select only subset of genes to see on plot
  if(!is.null(corr_genes_selection)) correlations <- correlations[rownames(correlations) %in% corr_genes_selection,]
  
  arrows = data.frame(x1 = rep(0, nrow(correlations)), y1 = rep(0, nrow(correlations)), x2 = correlations[, firstPC], y2 = correlations[,secondPC])
  
  ## PCA plot ##
  scores$clust <- factor(scores$clust)
  myrownames <- rownames(scores)
  
  pcaplot <- ggplot()
  pcaplot <- pcaplot + geom_point(data = scores, aes_string(x = firstPC, 
                                                            y = secondPC, colour = "clust"), size = 5)
  pcaplot <- pcaplot + scale_color_manual(values = colorClust)
  if (segments) {
    pcaplot <- pcaplot + geom_segment(data = scores, aes_string(x = firstPC, y = secondPC, xend = "xend", yend = "yend", colour = "clust"), 
                                      alpha = 0.3, size = 1)
  }
  if(centerpoint) {
    pcaplot <- pcaplot + geom_point(data = scores, aes(x=xend, y=yend), colour = "grey", alpha = 0.8, size = 16)
    pcaplot <- pcaplot + geom_point(data = scores, aes(x=xend, y=yend, colour = clust), alpha = 0.8, size = 14)
  }
  
  pcaplot <- pcaplot + theme_bw()
  pcaplot <- pcaplot + geom_hline(aes(yintercept = 0), colour = "lightgrey", 
                                  alpha = 0.7) + geom_vline(aes(xintercept = 0), colour = "lightgrey", 
                                                            alpha = 0.7)
  if (labeling) {
    pcaplot <- pcaplot + geom_text(data = scores, aes_string(x = firstPC, 
                                                             y = secondPC), label = myrownames, size = 5, colour = "darkgrey")
  }
  pcaplot <- pcaplot + xlab(bquote(atop(.(Xtitle), atop(italic(.(Xsubtitle)),""))))
  pcaplot <- pcaplot + ylab(bquote(atop(.(Ytitle), atop(italic(.(Ysubtitle)),""))))
  
  if (corcirc) {
    pcaplot <- pcaplot + geom_path(data = corcir, aes(x = x, y = y), colour = "gray65")
  }
  
  if (corGenes) {
    pcaplot <- pcaplot + geom_segment(data = arrows, aes(x = x1, y = y1, xend = x2, yend = y2), colour = "gray65")
    pcaplot <- pcaplot + geom_text(data = correlations, aes_string(x = firstPC,y = secondPC), label = rownames(correlations), size = 8)
  }
  
  pcaplot <- pcaplot + theme(legend.position = "none")
  pcaplot <- arrangeGrob(pcaplot)
  
  
  return(list(PCAdata = PCAdata, screeplot = screeplot, PCAscores = scores, pcaplot = pcaplot))
  
}


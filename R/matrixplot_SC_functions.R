
#' plotMarginal_SC function
#'
#' Function creates a histogram and density line for assay for matrixplot_SC()
#' @param variable vector of assay expression values
#' @param lwd linewidth
#' @return returns a histogram with line density
#' @export
#' @details NA
#' @references https://hlplab.wordpress.com/2012/03/20/correlation-plot-matrices-using-the-ellipse-library/
#' @examples
#' plotMarginal_SC()


plotMarginal_SC <- function(variable, lwd = 2) {
  
  variable <- variable[!is.na(variable)]
  density <- density(variable)
  h <- hist(variable, plot = FALSE)
  jitVar <- jitter(variable)
  yhigh <- max(max(h$density), max(density$y))
  ylow <- 0
  xticks <- pretty(c(variable, jitVar), min.n = 3)
  h <- hist(variable, freq = F, main = "", ylim = c(ylow, yhigh), 
            xlim = c(min(xticks), max(xticks)), xlab = "", ylab = " ", 
            axes = F, col = "grey", add = F)
  ax1 <- axis(1, line = 0.3, at = xticks, lab = xticks)
  par(las = 0)
  ax2 <- axis(2, at = c(0, max(max(h$density), max(density$y))/2, 
                        max(max(h$density), max(density$y))),
              labels = c("","Density", ""), lwd.ticks = 0, pos = NA, mgp = c(3, 0.2, 0), cex.axis = 1.7, mgp = c(3, 0.7, 0))
  rug(jitVar)
  lines(density$x[density$x >= min(ax1) & density$x <= max(ax1)], 
        density$y[density$x >= min(ax1) & density$x <= max(ax1)], 
        lwd = lwd)
  
}



#' plotScatter_SC function
#'
#' Function creates a scatterplot for matrixplot_SC()
#' @param dfr dataframe with expression values for pairwise comparison
#' @return returns a scatterplot between 2 variables
#' @export
#' @details NA
#' @references https://hlplab.wordpress.com/2012/03/20/correlation-plot-matrices-using-the-ellipse-library/
#' @examples
#' plotScatter_SC()
plotScatter_SC <- function(dfr, xVar, yVar, cexPoints = 1.3, cexXAxis = 1.3, cexYAxis = 1.3, lwd = 2) {
  
  compl_dfr <- dfr[complete.cases(dfr[, c(xVar, yVar)]), ]
  xVar <- compl_dfr[, xVar]
  yVar <- compl_dfr[, yVar]
  xticks <- pretty(xVar)
  yticks <- pretty(yVar)
  xlow <- min((min(xVar) - 0.1 * min(xVar)), min(xticks))
  xhigh <- max((max(xVar) + 0.1 * max(xVar)), max(xticks))
  ylow <- min((min(yVar) - 0.1 * min(yVar)), min(yticks))
  yhigh <- max((max(yVar) + 0.1 * max(yVar)), max(yticks))
  plot(xVar, yVar, col = "black", pch = 21, bg = "grey", ylab = "", 
       xlab = "", axes = F, ylim = c(ylow, yhigh), xlim = c(xlow,xhigh), cex = cexPoints)
  par(las = 1)
  axis(1, line = 0.4, labels = xticks, at = xticks, cex.axis = cexXAxis)
  axis(2, line = 0.2, labels = yticks, at = yticks, cex.axis = cexYAxis)
  fit <- vector("list", 4)
  fit[[1]] <- lm(yVar ~ poly(xVar, 1, raw = TRUE))
  fit[[2]] <- lm(yVar ~ poly(xVar, 2, raw = TRUE))
  fit[[3]] <- lm(yVar ~ poly(xVar, 3, raw = TRUE))
  fit[[4]] <- lm(yVar ~ poly(xVar, 4, raw = TRUE))
  Bic <- vector("numeric", 4)
  for (i in 1:4) {
    Bic[i] <- BIC(fit[[i]])
  }
  bestModel <- which.min(Bic)
  
  poly.pred <- function(fit) {
    f <- vector("character", 0)
    for (i in seq_along(coef(fit))) {
      if (i == 1) {
        temp <- paste(coef(fit)[[i]])
        f <- paste(f, temp, sep = "")
      }
      if (i > 1) {
        temp <- paste("(", coef(fit)[[i]], ")*", "x^", 
                      i - 1, sep = "")
        f <- paste(f, temp, sep = "+")
      }
    }
    x <- seq(min(xVar), max(xVar), length.out = 100)
    predY <- eval(parse(text = f))
    lines(x, predY, lwd = lwd)
  }
  
  poly.pred(fit[[bestModel]])
  
}



#' plotCorValue_SC()
#'
#' Function calculates pairwise correlation between 2 variables
#' @param xVar, yVar variables for correlation score
#' @return returns a spearman correlation score
#' @export
#' @details NA
#' @references https://hlplab.wordpress.com/2012/03/20/correlation-plot-matrices-using-the-ellipse-library/
#' @examples
#' plotCorValue_SC()
plotCorValue_SC <- function(xVar, yVar, cexText = 7.5, cexCI = 1.8, cor_use = "complete.obs", 
                            cor_method = "spearman") {
  
  plot(1, 1, type = "n", axes = FALSE, ylab = "", xlab = "")
  mycor <- function(x, y, my_use = cor_use, my_method = cor_method) {
    cor(x, y, use = my_use, method = my_method)
  }
  if (mycor(xVar, yVar) >= 0 & mycor(xVar, yVar) < 1) {
    lab = bquote(italic(r) == .(substr(x = formatC(round(mycor(xVar, 
                                                               yVar), 2), format = "f", digits = 2), start = 2, 
                                       stop = 4)))
    text(1, 1, labels = lab, cex = abs(mycor(xVar, yVar)) * 
           cexText, col = "#91cf60")
  }
  if (mycor(xVar, yVar) < 0) {
    lab = bquote(italic(r) == -.(substr(x = formatC(round(mycor(xVar, 
                                                                yVar), 2), format = "f", digits = 2), start = 3, 
                                        stop = 5)))
    text(1, 1, labels = lab, cex = abs(mycor(xVar, yVar)) * 
           cexText, col = "#fc8d59")
  }
  if (mycor(xVar, yVar) == 1) {
    lab = bquote(italic(r) == 1)
    text(1, 1, labels = lab, cex = abs(mycor(xVar, yVar)) * 
           cexText, col = "#1a9850")
  }
  ctest <- cor.test(xVar, yVar)
  CIlow <- formatC(round(ctest$conf.int[1], 2), format = "f", 
                   digits = 2)
  CIhigh <- formatC(round(ctest$conf.int[2], 2), format = "f", 
                    digits = 2)
  if (CIlow < 0) {
    CIlow <- paste("-", substr(CIlow, 3, 5), sep = "")
  }
  if (CIlow > 0) {
    CIlow <- substr(CIlow, 2, 4)
  }
  if (CIhigh < 0) {
    CIhigh <- paste("-", substr(CIhigh, 3, 5), sep = "")
  }
  if (CIhigh > 0) {
    CIhigh <- substr(CIhigh, 2, 4)
  }
  text(1, 0.8, labels = paste("95% CI: [", CIlow, ", ", CIhigh, 
                              "]", sep = ""), cex = cexCI)
  
}




#' matrixplot_SC()
#'
#' Function creates a pairwise matrixplot between assays and shows correlation, distribution and scatterplot
#' @param fluidSCproc fluidSCproc S3 object
#' @param based_on_values values to use, defaults to "log2Ex"
#' @param sorted_selected_genes subset of genes to select for visualization
#' @param not_detected_value value to be considered not expressed
#' @param replace_not_detected_with_NA boolean to replace not detected values with NA, so that they can be excluded at correlation calculation step
#' @param cor_method correlation algorithm method
#' @param use_method data to use for correlation algorithm
#' @return returns a histogram with line density
#' @export
#' @details NA
#' @references https://hlplab.wordpress.com/2012/03/20/correlation-plot-matrices-using-the-ellipse-library/
#' @examples
#' matrixplot_SC()
#' 
matrixplot_SC <- function(fluidSCproc,  based_on_values = "log2Ex",sorted_selected_genes = c("Pou5f1","Psma3","Trp53","Nanog","Cdx2"),
                          not_detected_value = 0, replace_not_detected_with_NA = F, cor_method = "spearman", use_method = "pairwise.complete.obs") {
  
  if(nargs() == 0) stop(paste0("you need to provide parameters, for more info see ?",sys.call()))
  
  normFluidCt <- fluidSCproc$data
  dataset <- dcast(normFluidCt, formula = Samples ~ Assays, value.var = based_on_values); rownames(dataset) <- dataset$Samples; dataset <- dataset[, -1]
  
  # replace not detected data with NA, can be excluded with selecting proper use method with cor function
  if(replace_not_detected_with_NA) dataset[dataset == not_detected_value] <- NA
  
  
  dataset <- dataset[ ,colnames(dataset) %in% sorted_selected_genes]
  dataset <- dataset[,sorted_selected_genes]
  
  variables <- colnames(dataset)
  l <- length(variables)
  if (l > 1) {
    par(mfrow = c(l, l), cex.axis = 1.3, mar = c(3, 4, 2, 1.5) + 0.1, oma = c(0, 2.2, 2, 0))
  }
  for (row in seq_len(l)) {
    for (col in seq_len(l)) {
      if (row == col) {
        plotMarginal_SC(dataset[[variables[row]]])
      }
      if (col > row) {
        plotScatter_SC(dataset, variables[col], variables[row])
      }
      if (col < row) {
        if (l < 7) {
          plotCorValue_SC(dataset[[variables[col]]], dataset[[variables[row]]], cexCI = 1.2, cor_use = use_method, cor_method = cor_method)
        }
        if (l >= 7) {
          plotCorValue_SC(dataset[[variables[col]]], dataset[[variables[row]]], cexCI = 1.2, cor_use = use_method, cor_method = cor_method)
        }
      }
    }
  }
  if (l == 1) {
    mtext(text = variables[1], side = 1, cex = 1.9, line = 3)
  }
  if (l > 1) {
    textpos <- seq(1/(l * 2), (l * 2 - 1)/(l * 2), 2/(l * 
                                                        2))
    for (t in seq_along(textpos)) {
      mtext(text = variables[t], side = 3, outer = TRUE, at = textpos[t], cex = 1.9, line = -0.8)
      mtext(text = variables[t], side = 2, outer = TRUE, at = rev(textpos)[t], cex = 1.9, line = -0.1)
    }
  }
  
}

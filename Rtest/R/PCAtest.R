#' LEMC: A likelihood-based method for estimating methylation-induced correlation among human complex traits
#' @title PCA.analysis
#'
#' @description This function is for PCA analysis to transform the methylation data.
#'
#' @details you can use this function to calculate the eigen value and eigen matrix of data.
#'
#' @param data.omics data.omics is input data frame or matrix
#'
#' @import Morpho
#'
#' @return a dataframe
#'
#' @examples
#' \dontrun{
#' ## a simulated DNA methylation data
#' load("data/data.example.Rdata")
#' library(Morpho)
#' PCA.data <- PCA.analysis(data.omics)
#' PCA.data
#' }
#' @export

### PCA of omics data
PCA.analysis <- function(data.omics, lambda = 1) {
  res.pca <- prcompfast(data.omics, retx = TRUE, center = TRUE, scale. = TRUE)
  omics.pca <- res.pca$x

  comprob <- summary(res.pca)
  omics.pca.top <- omics.pca[, 1:sum(comprob$importance[1, ] > lambda)]
  X <- data.frame(scale(omics.pca.top))
  return(X)
}

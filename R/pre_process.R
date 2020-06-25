#' @title Pre-processing gene expression matrix
#'
#' @description Use this function to perform a number of pre-processing steps, including
#' log transformation, normalization and standardization of gene expression matrix.
#' NOTE: If multiple options are set true, the sequence of pre-processing steps
#' performed are variance filtering, log2 transformation, upper quantile normalization
#' and standardization.
#'
#' @param mat a gene expression matrix with genes in columns and samples in rows
#' @param var.filter logical. specify whether genes are filtered based on variance. Default = TRUE
#' @param log logical. specify whether log2(x+1) transformation should be performed. Default = FALSE
#' @param std logical. specifiy whether standardization should be performed. Default = TRUE
#' @param UQ logical. specify whether upper quantile normalization should be performed. Default = FALSE
#'
#' @return a processed gene expression matrix
#'
#' @examples
#' Proc.GA.pcg <- pre_process(GA.pcg)
#'
pre_process <- function (mat, var.filter = TRUE, log = FALSE, UQ = FALSE, std = FALSE) {
  if (mode(mat)!="numeric" | !is(mat, "matrix")) stop ("Error: input data must be a numeric matrix")

  if (var.filter == TRUE) {
    gvar <- apply(mat, 2, var)
    mat <- mat[, which(gvar!=0)]
  }

  if(log == TRUE) mat <- log2(mat + 1)

  if (UQ == TRUE) mat <- apply(mat, 2, function (x) x/quantile(x, 0.75))

  if (std == TRUE) mat <- scale(mat)

  return(mat)
}

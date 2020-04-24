#' @title Select informative PCG subset
#'
#' @description Function to filter informative protein coding genes based on correlation with circRNA expression
#'
#' @param train.pcg training protein coding dataset,
#' which should be a numeric matrix with row names indicating samples,
#' and column names indicating protein coding gene IDs.
#' @param train.circ training circRNA expression dataset,
#' which should be a numeric matrix with row names indicating samples,
#' and column names indicating circRNA IDs.
#' @param gene.index either gene name (character) or index (column number) of circRNA to be imputed.
#' @param num number of informative protein coding genes to be used in constructing imputation model. Default is 50 genes.
#'
#' @return a numeric matrix, which should be a subset of protein coding genes correlated with circRNA of interest.
corf <- function (train.pcg, train.circ, gene.index, num = 50) {
  pcor <- abs(cor(train.pcg, train.circ[, gene.index]))
  r_pcor <- rank(-pcor)
  gin <- which (r_pcor < num)
  temp_pcg <- train.pcg[, gin]
  return(temp_pcg)
}

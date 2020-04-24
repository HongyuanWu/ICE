#' @title Match gene expression matrices by gene ID
#'
#' @description Returns expression matrices with common genes
#'
#' @param train.pcg Gene expression matric with gene IDs as column names
#' @param new.pcg Gene expression matric with gene IDs as column names
#'
#' @return A list containing two matrices with common genes
#'
#' @examples
#' temp <- match_mat(GA.pcg, HS.pcg)
#' Matched.GA.pcg <- data.matrix(temp[[1]])
#' Matched.HS.pcg <- data.matrix(temp[[2]])
#'
#' @export
match_mat <- function (train.pcg, new.pcg) {

  if (mode(train.pcg) != "numeric" | class(train.pcg) != "matrix" |
      mode(new.pcg) != "numeric" | class(new.pcg) != "matrix" ) stop ("Error: input data must be a numeric matrix")

  y <- match (colnames(train.pcg), colnames(new.pcg))
  y1 <- which(!is.na(y))
  y2 <- na.omit(y)

  ncom <- length(y1)
  if (ncom < 0.25 * nrow(train.pcg)) stop ("Error: < 25% common genes in training and test dataset")
  if (ncom < 200) stop ("Error: < 200 common genes in the two datasets")

  new.train.pcg <- train.pcg[, y1]
  new.new.pcg <- new.pcg[, y2]
  cmat <- list(new.train.pcg, new.new.pcg)

  return(cmat)
}

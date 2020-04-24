#' ICE: Impute CircRNA Expression
#'
#' @docType package
#' @name ICE
#'
#' @import randomForest
#' @import FNN
#' @import e1071
NULL

#' ICE: Impute expression of a certain series of circRNAs of interest
#'
#' @param train.circ training circRNA expression dataset, which should be a numeric matrix with row names indicating
#' samples, and column names indicating circRNA IDs.
#' @param train.pcg training protein coding dataset, which should be a numeric matrix with with row names indicating
#' samples, and column names indicating protein coding gene IDs.
#' @param new.pcg protein coding expression dataset to be used for prediction, which should be a numeric matrix with row names indicating
#' samples, and column names indicating protein coding gene IDs.
#' @param gene.index either gene name (character) or index (column number) of circRNA to be imputed.
#' @param method method for imputation, either "RF" for random forests, "KNN" for K-nearest neighbor or
#' "SVM" for support vector machines. Uses KNN by default.
#' @param num number of informative protein coding genes to be used in constructing imputation model.
#' Default is 50 genes.
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' data(iMIRAGE.datasets)
#' imirage(GA.pcg, GA.mir, HS.pcg, gene_index="hsa-let-7c", method="KNN", num=50)
#' imirage(GA.pcg, GA.mir, HS.pcg, gene_index=25, method="KNN", num=50)
ICE <- function(train.circ, train.pcg, new.pcg, gene.index, method = "KNN", num = 50, ...) {
  if (mode(train.pcg) != "numeric" | mode(train.circ) != "numeric" | mode(new.pcg) != "numeric" |
      class(train.pcg) != "matrix" | class(train.circ) != "matrix" | class(new.pcg) != "matrix") stop ("Error: input data must be a numeric matrix")


  if (mode(gene.index) == "numeric" & gene.index > ncol(train.circ)) stop ("Error: miRNA not found in training dataset")
  if (mode(gene.index) == "character" & is.na (match (gene.index, colnames(train.circ)))) stop ("Error: miRNA not found")

  temp <- match_mat(train.pcg, new.pcg)
  train.pcg <- temp[[1]]
  new.pcg <- temp[[2]]

  train.pcg <- pre_process(train.pcg)
  train.circ <- filter_circ(train.circ)

  if (mode(gene.index) == "numeric") y <- train.circ[, gene.index]
  if (mode(gene.index) == "character") y <- train.circ[, match(gene.index, colnames(train.circ))]

  if (sd(y)==0) stop ("Error: Standard deviation of circRNA is 0")
  if (nrow(train.pcg) < 50) warning("Warning: Sample size is < 50")

  if (num >= 50 & ncol(train.pcg) >=50 ) x <- corf(train.pcg, train.circ, gene.index, num)
  if (ncol(train.pcg) < 50) x <- train.pcg

  if (method == "RF") {
    rfit <- randomForest(x, y, ntree=250, ...)
    predict.y <-  predict(rfit, new.pcg)
    return(predict.y)
  }

  if (method == "KNN") {
    mX <- match(colnames(x), colnames(new.pcg))
    predict.y <- knn.reg(train = x, test = new.pcg[,mX], y = y, k=50)$pred
    return(predict.y)
  }

  if (method == "SVM") {
    mX <- match(colnames(x), colnames(new.pcg))
    svm.fit <- svm(x = x, y = y)
    predict.y <- predict(svm.fit, new.pcg[,mX])
    return(predict.y)
  }
}

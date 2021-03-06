#' ICE: Impute CircRNA Expression
#'
#' @docType package
#' @name ICE
#'
#' @import randomForest
#' @import FNN
#' @import e1071
#' @importFrom glmnet glmnet predict.glmnet
#' @importFrom caret train
#' @importFrom pls pcr
#' @importFrom snowfall sfInit sfGetCluster sfExport sfStop
#' @importFrom pbapply pblapply pboptions
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
#' Default is 100 genes.
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' data("mionco.circ")
#' data("mionco.pcg")
#' pred.circ <- ICE(train.pcg = mionco.pcg, train.circ = mionco.circ, new.pcg = ccle.pcg, gene.index = "hsa_circ_0000801")
ICE <- function(train.pcg, train.circ, new.pcg, gene.index, method = "KNN", num = 50, filter = TRUE, ...) {
  if (mode(train.pcg) != "numeric" | mode(train.circ) != "numeric" | mode(new.pcg) != "numeric" |
      !is(train.pcg, "matrix") | !is(train.circ, "matrix") | !is(new.pcg, "matrix")) stop ("Error: input data must be a numeric matrix")


  if (mode(gene.index) == "numeric" & gene.index > ncol(train.circ)) stop ("Error: circRNA not found in training dataset")
  if (mode(gene.index) == "character" & is.na (match (gene.index, colnames(train.circ)))) stop ("Error: circRNA not found")

  temp <- match_mat(train.pcg, new.pcg)
  train.pcg <- temp[[1]]
  new.pcg <- temp[[2]]

  if(filter) {
    train.pcg <- pre_process(train.pcg)
    train.circ <- scale(filter_circ(train.circ))
  }

  if (mode(gene.index) == "numeric") y <- train.circ[, gene.index]
  if (mode(gene.index) == "character") y <- train.circ[, match(gene.index, colnames(train.circ))]

  if (sd(y)==0) stop ("Error: Standard deviation of circRNA is 0")
  if (nrow(train.pcg) < num) warning("Warning: Sample size is less than desired feature number")

  if (ncol(train.pcg) >= num & method != 'lasso' & method != "EN") x <- corf(train.pcg, train.circ, gene.index, num)
  if (ncol(train.pcg) < num) x <- train.pcg

  if (method == "RF") {
    rfit <- randomForest(x, y, ntree=500, mtry=33, nodesize=1)
    predict.y <- predict(rfit, new.pcg)
    return(predict.y)
  }

  if (method == "KNN") {
    mX <- match(colnames(x), colnames(new.pcg))
    predict.y <- knn.reg(train = x, test = new.pcg[,mX], y = y, k = 1)$pred
    return(predict.y)
  }

  if (method == "SVM") {
    mX <- match(colnames(x), colnames(new.pcg))
    svm.fit <- svm(x = x, y = y, cost = 2222, gamma = 0.12, type = 'eps-regression', kernel = 'radial')
    predict.y <- predict(svm.fit, new.pcg[,mX])
    return(predict.y)
  }

  if (method == "lasso") {
    # https://stackoverflow.com/a/32185491
    lasso.model <- glmnet(x, y, alpha = 1, family = 'gaussian')
    predict.y <- predict.glmnet(lasso.model, s = 0.01530921, newx = x)
    return(predict.y)
  }

  if (method == "EN") {
    # http://www.sthda.com/english/articles/37-model-selection-essentials-in-r/153-penalized-regression-essentials-ridge-lasso-elastic-net/
    en.model <- glmnet(x, y, alpha = 0.9909091, family = 'gaussian')
    predict.y <- predict.glmnet(lasso.model, s = 0.03019828, newx = x)
    return(predict.y)
  }

  if (method == 'PCR') {
    # http://www.milanor.net/blog/performing-principal-components-regression-pcr-in-r/
    temp.df <- data.frame(x, y)
    pcr.model <- pcr(y ~ ., data = temp.df, scale =TRUE, validation = "CV")
    predict.y <- predict(pcr.model, x, ...)[, , ]
    return(predict.y)
  }
}

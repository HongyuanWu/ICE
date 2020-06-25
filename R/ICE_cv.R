#' @title Cross validation function for ICE imputation accuracy
#'
#' @description Function to obtain accuracy parameters: correlation coefficient, P-value and RMSE of
#' imputation model
#'
#' @param train.circ training circRNA expression dataset, which should be a numeric matrix with row names indicating
#' samples, and column names indicating circRNA IDs.
#' @param train.pcg training protein coding dataset, which should be a numeric matrix with with row names indicating
#' samples, and column names indicating protein coding gene IDs.
#' @param gene.index either gene name (character) or index (column number) of circRNA to be imputed.
#' @param method method for imputation, either "RF" for random forests, "KNN" for K-nearest neighbor or
#' "SVM" for support vector machines.
#' @param num number of informative protein coding genes to be used in constructing imputation model.
#' Default is 100 genes.
#' @param folds number specifying folds (k) of cross validation to obtain imputation accuracy.
#' Default is k=10.
#' @param ... optional parameters that can be passed on to the machine-learning functions
#' RF (\link[randomForest]{randomForest}), KNN (\link[FNN]{knn.reg}) or SVM(\link[e1071]{svm})
#' @param filter if training datasets should be filtered? Default is TRUE, but should be set to FALSE when called in loop.
#' @param cor.method a character string indicating
#' which correlation coefficient is to be used for the correlation test.
#' Default is "spearman", could also be "pearson".
#'
#' @return a matrix with three values corresponding to correlation coefficient,
#' P-value of the fit and root mean squared error (RMSE).
#'
#' @examples
#' data(iMIRAGE.datasets)
#' ICE_cv(GA.pcg, GA.mir, gene.index="hsa-let-7c", method="KNN", num=100)
#' ICE_cv(GA.pcg, GA.mir, gene.index=25, method="KNN", num=100)
#'
#' @export ICE_cv
ICE_cv <- function (train.pcg, train.circ, gene.index, num = 50,
                    method = "KNN", folds = 10, filter = TRUE, cor.method = 'spearman', ...) {

  if (mode(gene.index) != "numeric" & mode(gene.index) != "character") stop ("Error: circRNA not found in training dataset. Please check gene name or rownumber")
  if (mode(gene.index) == "numeric" & gene.index > ncol(train.circ))  stop ("Error: circRNA not found in training dataset. Please check ID or rownumber")
  if (mode(gene.index) == "character" & is.na (match (gene.index, colnames(train.circ)))) stop ("Error: circRNA not found. Please check ID or rownumber")

  if (filter) {
    train.pcg <- pre_process(train.pcg)
    train.circ <- scale(filter_circ(train.circ))
  }

  if (ncol(train.pcg) >= num & method != 'lasso' & method != "EN") train.pcg <- scale(corf(train.pcg, train.circ, gene.index, num))
  if (ncol(train.pcg) < num) train.pcg <- scale(train.pcg)

  cv.res <- matrix(nrow=folds, ncol=3)
  colnames (cv.res) <- c("CC", "P-Value", "RMSE")

  RS <- nrow(train.pcg)
  while (RS %% folds != 0) {
    RS = RS - 1
  }
  kgrp <- split(sample(1:RS, RS, replace=F), 1:folds)

  for (k in 1:folds) {

    ind <- unlist(kgrp[[k]])
    x <- train.pcg [-ind, ]
    testx <- train.pcg [ind, ]

    if (mode(gene.index)=="numeric") {
      y <- train.circ[-ind, gene.index]
      actual.y <- train.circ [ind, gene.index]
    }

    if (mode(gene.index)=="character") {
      y <- train.circ[-ind, match(gene.index, colnames(train.circ))]
      actual.y <- train.circ[ind, match(gene.index, colnames(train.circ))]
    }

    if(method=="RF") {
      #randomForest
      imp.rf <- randomForest(x, y, ntree=500, mtry=33, nodesize=1)
      predict.y <- predict(imp.rf, testx)
      r.rf <- suppressWarnings(cor.test(predict.y, actual.y, method = cor.method))
      rmse.rf <- sqrt(mean((actual.y-predict.y)^2))
      cv.res[k, 1:3] <- c(r.rf$estimate, r.rf$p.value, rmse.rf)
      next()
    }

    if (method=="KNN") {
      #knn regression
      knn.fit <- knn.reg(train = x, test = testx, y = y, k = 1, ...)
      r.knn <- suppressWarnings(cor.test(knn.fit$pred, actual.y, method = cor.method))
      rmse.knn <- sqrt(mean((actual.y - knn.fit$pred)^2))
      cv.res[k, 1:3] <- c(r.knn$estimate, r.knn$p.value, rmse.knn)
      next()
    }

    if (method=="SVM") {
      #svm regression
      imp.svm <- svm(x, y, cost = 2222, gamma = 0.12, type = 'eps-regression', kernel = 'radial')
      predict.y <- predict(imp.svm, testx)
      r.svm <- suppressWarnings(cor.test(predict.y, actual.y, method = cor.method))
      rmse.svm <- sqrt(mean((actual.y-predict.y)^2))
      cv.res[k, 1:3] <- c(r.svm$estimate, r.svm$p.value, rmse.svm)
      next()
    }

    if (method=="lasso") {
      lasso.model <- glmnet(x, y, alpha = 1, family = 'gaussian')
      predict.y <- predict.glmnet(lasso.model, s = 0.01530921, newx = testx)
      r.lasso <- suppressWarnings(cor.test(predict.y, actual.y, method = cor.method))
      rmse.lasso <- sqrt(mean((actual.y - predict.y) ^ 2))
      cv.res[k, 1:3] <- c(r.lasso$estimate, r.lasso$p.value, rmse.lasso)
      next()
    }

    if (method == "EN") {
      en.model <- glmnet(x, y, alpha = 0.9909091, family = 'gaussian')
      predict.y <- predict.glmnet(en.model, s = 0.03019828, newx = testx)
      r.en <- suppressWarnings(cor.test(predict.y, actual.y, method = cor.method))
      rmse.en <- sqrt(mean((actual.y - predict.y) ^ 2))
      cv.res[k, 1:3] <- c(r.en$estimate, r.en$p.value, rmse.en)
      next()
    }

    if (method == 'PCR') {
      temp.df <- data.frame(x, y)
      pcr.model <- pcr(y ~ ., data = temp.df, scale =TRUE, validation = "CV", ...)
      predict.y <- predict(pcr.model, testx, ...)[, , ]
      r.pcr <- suppressWarnings(cor.test(predict.y, actual.y, method = cor.method))
      rmse.pcr <- sqrt(mean((actual.y - predict.y) ^ 2))
      cv.res[k, 1:3] <- c(r.pcr$estimate, r.pcr$p.value, rmse.pcr)
      next()
    }
  }

  return (cv.res)
}

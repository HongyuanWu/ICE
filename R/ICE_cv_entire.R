#' @title ICE cross-validation loop function for full circRNA matrix
#' @description Convinient wrapper for \link[ICE]{ICE_cv} that performs cross-validation analysis for
#' assessing imputation accuracies for all circRNAs using the training datasets
#' @param train.pcg training protein coding dataset. a numeric matrix with with row names indicating
#' samples, and column names indicating protein coding gene IDs.
#' @param train.circ training circRNA expression dataset. a numeric matrix with row names indicating
#' samples, and column names indicating circRNA IDs.
#' @param method method for imputation, either "RF" for random forests, "KNN" for K-nearest neighbor or
#' "SVM" for support vector machines.
#' @param ... optional parameters that can be passed on to the machine-learning method:
#' RF (\link[randomForest]{randomForest}), KNN (\link[FNN]{knn.reg}) or SVM(\link[e1071]{svm})
#'
#' @examples
#' ICE_cv_entire(train.pcg = GA.pcg, train.circ = GA.mir, method = "KNN")
#' ICE_cv_entire(train.pcg = GA.pcg, train.circ = GA.mir, method = "SVM")
#'
#' @return a matrix containing Spearman's correlation coefficient, P-value and RMSE from the cross-validation analysis
#' of the complete circRNA training dataset
#' @export ICE_cv_entire
ICE_cv_entire <- function (train.pcg, train.circ, method = "KNN", ...) {
  cv.loop <- list()
  for (i in 1:ncol(train.circ)) {
    cv.loop[[i]] <- ICE_cv(train.pcg, train.circ, gene.index = i, method = method, filter = FALSE, ...)
    print(i)
  }

  cv.results <- cv_proc(cv.loop)
  return(cv.results)
}

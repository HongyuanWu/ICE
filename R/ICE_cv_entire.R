#' @title ICE cross-validation for entire circRNA matrix
#' @description Convinient wrapper for \link[ICE]{ICE_cv} that performs cross-validation analysis for
#' assessing imputation accuracies for all circRNAs using the training datasets
#'
#' @param train.pcg training protein coding dataset. a numeric matrix with with row names indicating
#' samples, and column names indicating protein coding gene IDs.
#' @param train.circ training circRNA expression dataset. a numeric matrix with row names indicating
#' samples, and column names indicating circRNA IDs.
#' @param method method for imputation, either "RF" for random forests, "KNN" for K-nearest neighbor or
#' "SVM" for support vector machines.
#' @param ... optional parameters that can be passed on to the machine-learning method:
#' RF (\link[randomForest]{randomForest}), KNN (\link[FNN]{knn.reg}) or SVM(\link[e1071]{svm})
#' @param cor.method a character string indicating
#' which correlation coefficient is to be used for the correlation test.
#' Default is "spearman", could also be "pearson".
#'
#' @examples
#' ICE_cv_entire(train.pcg = GA.pcg, train.circ = GA.mir, method = "KNN")
#' ICE_cv_entire(train.pcg = GA.pcg, train.circ = GA.mir, method = "SVM")
#'
#' @return a matrix containing Spearman's correlation coefficient, P-value and RMSE from the cross-validation analysis
#' of the complete circRNA training dataset
#' @export ICE_cv_entire
ICE_cv_entire <- function (train.pcg, train.circ, method = "KNN", cor.method = 'spearman', cutoff = 75, folds = 10, ncores = 16, filter = TRUE, ...) {

  if (method != 'KNN' & method != 'PCR' & ncores != 1) {
    sfInit(parallel = TRUE, cpu = ncores)
    cl <- sfGetCluster()
  }


  if (filter) {
    train.pcg <- pre_process(train.pcg)
    # Here it is recommended to use parameters "cutoff = 95, threshold = 0"
    # for there may be cases that all the expression value of a certain circRNA
    # among a series samples used in one fold is equal to zero, which will bring Nans after scale() during pcr().
    train.circ <- scale(filter_circ(train.circ, cutoff = cutoff))
  } else {
    train.circ <- scale(train.circ)
  }


  if (method != 'KNN' & method != 'PCR' & ncores != 1) {
    sfExport('train.pcg', 'train.circ')
  }

  if (method != 'KNN' & method != 'PCR' & ncores != 1) {
    opb <- pboptions(use_lb = TRUE)
    on.exit(pboptions(opb))
    entire.cv <- pblapply(seq_len(NCOL(train.circ)),
                          function(x) ICE_cv(train.pcg, train.circ, gene.index = x,
                                             method = method, filter = FALSE, folds = folds, cor.method = cor.method, ... = ...), cl = cl)
  }else{
    entire.cv <- pblapply(seq_len(NCOL(train.circ)),
                          function(x) ICE_cv(train.pcg, train.circ, gene.index = x,
                                             method = method, filter = FALSE, folds = folds, cor.method = cor.method, ... = ...))
  }

  if (method != 'KNN' & method != 'PCR' & ncores != 1) {
    sfStop()
  }

  (cv.res <- t(vapply(entire.cv, function(x) apply(x, 2, mean), numeric(3))))
}

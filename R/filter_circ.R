#' @title Pre-standardization filter to select circRNA expressed in a specific proportion of samples
#'
#' @description return a subset of the circRNA expression dataset containing circRNAs that are expressed
#' in at least a specified percentage of samples based on a user-defined expression threshold
#'
#' @param train.circ an expression matrix with circRNA in columns, samples in rows
#' @param cutoff percentage of samples in which the circRNA should be expressed
#' @param threshold the numeric threshold defining expression of circRNA. Default threshold = 0.
#' @return filtered circRNA expression matrix
#'
#' @examples
#' New.mir <-  filter_circ(GA.mir, cutoff=95, threshold=0)
#'
filter_circ <- function (train.circ, cutoff = 75, threshold = 0) {

  if (mode(train.circ)!="numeric" | class(train.circ)!="matrix") stop ("Error: input data must be a numeric matrix")

  index <- 0

  cut <- cutoff*nrow(train.circ)/100

  for (i in 1:ncol(train.circ)) {
    if (sum(train.circ[, i] > threshold) > cut)
      index <- append(index, i)
  }

  index <- index[2:length(index)]
  f_train.circ <- train.circ[, index]

  return (f_train.circ)
}

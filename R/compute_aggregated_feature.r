#' Compute patient/tumor-level aggregated cell features
#'
#' Compute the min, max, 25th, 50th, 75th, and stdev of each cellular features across all tumors/cases
#'
#' @param data A data.frame contains cells in rows and features in columns;
#' the column 'tumor_ids' represents the tumor identifiers
#' @param output_dir A file.path string of output directory
#' @return A data.frame contains tumor-level feature summaries saved in
#'   \code{output_dir}
#' @importFrom stats quantile sd
#' @import utils
#' @export
compute_aggregated_feature <- function (data=NULL, output_dir = getwd()){

  ## ======================
  ##  compute statistics grouping by tumor_ids
  ## ======================
  tumor_ids <- NULL
  Q1_f <- plyr::ddply(data, plyr::.(tumor_ids), plyr::numcolwise(quantile,probs=0.25,na.rm = TRUE))
  rownames(Q1_f) <- Q1_f$tumor_ids; Q1_f$tumor_ids <- NULL
  colnames(Q1_f) <- paste0('Q1.',colnames(Q1_f))

  Q2_f <- plyr::ddply(data, plyr::.(tumor_ids), plyr::numcolwise(quantile,probs=0.5,na.rm = TRUE))
  rownames(Q2_f) <- Q2_f$tumor_ids; Q2_f$tumor_ids <- NULL
  colnames(Q2_f) <- paste0('Q2.',colnames(Q2_f))

  Q3_f <- plyr::ddply(data, plyr::.(tumor_ids), plyr::numcolwise(quantile,probs=0.75,na.rm = TRUE))
  rownames(Q3_f) <- Q3_f$tumor_ids; Q3_f$tumor_ids <- NULL
  colnames(Q3_f) <- paste0('Q3.',colnames(Q3_f))

  min_f <- plyr::ddply(data, plyr::.(tumor_ids), plyr::numcolwise(min,na.rm = TRUE))
  rownames(min_f) <- min_f$tumor_ids; min_f$tumor_ids <- NULL
  colnames(min_f) <- paste0('min.',colnames(min_f))

  max_f <- plyr::ddply(data, plyr::.(tumor_ids), plyr::numcolwise(max,na.rm = TRUE))
  rownames(max_f) <- max_f$tumor_ids; max_f$tumor_ids <- NULL
  colnames(max_f) <- paste0('max.',colnames(max_f))

  stdev_f <- plyr::ddply(data, plyr::.(tumor_ids), plyr::numcolwise(sd,na.rm = TRUE))
  rownames(stdev_f) <- stdev_f$tumor_ids; stdev_f$tumor_ids <- NULL
  colnames(stdev_f) <- paste0('stdev.',colnames(stdev_f))


  ## ======================
  ##  combine summary statistics
  ## ======================
  summaries_data <- merge(Q1_f,Q2_f, by='row.names')
  rownames(summaries_data) <- summaries_data$Row.names; summaries_data$Row.names <- NULL
  summaries_data <- merge(summaries_data,Q3_f, by='row.names')
  rownames(summaries_data) <- summaries_data$Row.names; summaries_data$Row.names <- NULL
  summaries_data <- merge(summaries_data,min_f, by='row.names')
  rownames(summaries_data) <- summaries_data$Row.names; summaries_data$Row.names <- NULL
  summaries_data <- merge(summaries_data,max_f, by='row.names')
  rownames(summaries_data) <- summaries_data$Row.names; summaries_data$Row.names <- NULL
  summaries_data <- merge(summaries_data,stdev_f, by='row.names')
  rownames(summaries_data) <- summaries_data$Row.names; summaries_data$Row.names <- NULL
  #ncol(summaries_data)/6

  ## ======================
  ##  saving results
  ## ======================
  setwd(output_dir)
  write.csv(summaries_data, file = paste0('aggregate_patient_features.csv'),row.names = TRUE)

}

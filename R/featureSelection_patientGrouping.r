#' Feature selection - patient grouping and prognostic significance
#'
#' Feature selection - patient grouping and prognostic significance:
#' using cell features obtained at each 10% increment of detection rates
#' (determined by \code{select_prognostic_feature}), tumors are divided into
#' groups which are evaluated for prognostic significance.
#' @param cell_feature_data A data.frame containing tumor-average cellular feature values: tumors in rows, features in columns
#' @param feature_data A data.frame containing the feature selection of each feature (rows) in each Cox PH Lasso run (columns)
#' @param surv_data A data.frame: tumors in rows (tumor_ids in row.names), with two colums "time" and "cens"
#' @param log10_featureNms A vector of characters of feature names whose values to be log10-transformed prior to scaling (z-score)
#' @param neighbor_size_byFrac A fraction (0 to 1) indicating the fraction of total samples/tumors to be
#' used in PhenoGraph
#' @param plot_bnm A character string for output plot filename identifier
#' @param output_dir A file.path string of output directory
#' @param seed An integer of seed (for repeatability)
#' @param detection_rate_cutoff A vector of fractions (0 to 1) indicating the (prognostically selected)
#' feature detection rates, to be tested
#' @param min_cluster_size An integer indicating the minimum cluster size (i.e. no. of cases/tumors) for which
#' a tumor subtype containing tumors less than this number will be excluded in the survival analysis
#' in \code{surv_analysis_fn} and in the corresponding output plot \code{Viz_heatmap_bySubtype}
#' @import Rphenograph ggplot2 survival survminer
#' @importFrom grDevices dev.off pdf
#' @importFrom stats relevel terms
#' @import utils
#' @export
featureSelection_patientGrouping <- function(cell_feature_data=NULL, feature_data= NULL,surv_data= NULL,
                                   log10_featureNms=NULL,
                                   neighbor_size_byFrac=0.05,
                                   plot_bnm='', output_dir = getwd(), seed=999,
                                   detection_rate_cutoff=seq(0,0.9,by=0.1),min_cluster_size=20){

  ## ======================
  ## parameters documentation
  ## ======================
  setwd(output_dir)
  capture.output(file = 'FS_patientGrouping_log.txt',
                 cat('neighbor_size_byFrac=',neighbor_size_byFrac,
                     '\nseed=',seed,
                     '\nlog10_featureNms=',paste0(log10_featureNms,collapse = '\t')))
  ## ======================
  ## calculate % feature detection
  ## ======================
  ## remove features with 0 time selection
  rs <- rowSums(feature_data[, colnames(feature_data)!='terms']==0)
  #sum(rs== (ncol(feature_data)-1) )
  all_runs2 <- feature_data[rs< (ncol(feature_data)-1),]
  cat('no. of feature remained at least 1 time = ', nrow(all_runs2), '\n')

  rownames(all_runs2) <- all_runs2$terms
  all_runs2$terms <- NULL

  detection_rate <- apply(X=all_runs2, MARGIN = 1, FUN = function(z){
    sum(abs(z)>0)/length(z)
  })

  detection_rate_df <- data.frame(features=names(detection_rate), percent_detect=detection_rate)

  ## sort by highest rates on top
  detection_rate_df <- detection_rate_df[order(detection_rate_df$percent_detect, decreasing = TRUE), ]
  #head(detection_rate_df)
  write.csv(detection_rate_df, file = file.path(output_dir,'feature_detection_rate.csv'), row.names = FALSE)
  prev_feature <- c(detection_rate_df$features,0)
  prev_feature2 <- prev_feature
  for (detection_rate in detection_rate_cutoff){
    detection_rate_df_sub <- detection_rate_df[detection_rate_df$percent_detect > detection_rate, ,drop=FALSE]
    if(sum(detection_rate_df_sub$features %in% prev_feature) == length(prev_feature))next
    if(nrow(detection_rate_df_sub)==0){
      warning('@detection rate=', detection_rate, ': no feature remained\n')
      break
    }
    cat('detection_rate=',detection_rate, '\n')
    prev_feature <- as.character(detection_rate_df_sub$features)
    feature_nms <- sapply(strsplit(as.character(detection_rate_df_sub$features), split = 'poly\\('), "[[", 2)
    feature_nms <- sapply(strsplit(feature_nms, split = '\\,'), "[[", 1)
    feature_nms <- data.frame(features=unique(feature_nms))
    if( nrow(feature_nms) == length(prev_feature2)){prev_feature2 <- feature_nms$features; next}
    prev_feature2 <- feature_nms$features
    output_dir2 <- file.path(output_dir, paste0('detection_rate',detection_rate))
    dir.create(output_dir2)
    write.csv(feature_nms, file = file.path(output_dir2,'selectedAggFeatures.csv'), row.names = FALSE)
    ############################
    ## data subsetting for selected features
    ############################
    data_subset <- cell_feature_data[, colnames(cell_feature_data) %in% feature_nms$features, drop=FALSE]
    rs_NA <- rowSums(is.na(data_subset))
    #sum(rs_NA>0)
    #View(data_subset[rs_NA>0, ])
    data_subset <- data_subset[rs_NA==0, , drop=FALSE]


    ## ======================
    ## data transformation
    ## ======================
    for (cc in 1:ncol(data_subset)){
      feature_cc <- colnames(data_subset)[cc]
      if(feature_cc %in% log10_featureNms){
        cat('LOG10:', feature_cc, '\n')
        data_subset[,cc] <-log10_transform(data_subset[,cc], method = 'log10')
      }
      data_subset[,cc] <- scale(data_subset[,cc])
    }
    #summary(data_subset)

    ## ======================
    ## subsetting surv data
    ## ======================
    surv_data_sub <- surv_data[match(rownames(data_subset), rownames(surv_data)),]
    stopifnot(identical(rownames(data_subset), rownames(surv_data_sub)))

    ## ======================
    ## Phenograph clustering of tumors: subtyping
    ## ======================
    neighbor_size <- floor(neighbor_size_byFrac * nrow(data_subset))
    Rphenograph_out <- Rphenograph(data_subset, k = neighbor_size)
    missing_case <- setdiff(1:nrow(data_subset),names(igraph::membership(Rphenograph_out[[2]])))
    if(length(missing_case) >0)data_subset <- data_subset[-missing_case, ]

    tumor_clusters <- data.frame(tumor_ids= rownames(data_subset), cluster_no=Rphenograph_out[[2]]$membership)
    write.csv(tumor_clusters, file = file.path(output_dir2,'FS_patientGrouping.csv'), row.names = FALSE)

    ## ======================
    ## surv analysis using identified tumor subtypes
    ## ======================
    surv_analysis_fn(surv_data=surv_data, subtype=tumor_clusters,out_bnm='FS_patientGrouping',
                     output_dir=output_dir2,min_cluster_size= min_cluster_size, KM_xlab='Survivals')

    data_subset2 <- data_subset
    data_subset2$tumor_ids <- rownames(data_subset2)
    Viz_heatmap_bySubtype(data = data_subset2, subtype=tumor_clusters,
                             output_dir=output_dir2,plot_bnm='featureSelection',
                             min_cluster_size=min_cluster_size,
                             value_name='Value', row_title='Feature')
  }#end detection_rate

}

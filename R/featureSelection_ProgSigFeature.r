#' Prognostic feature selection using Cox PH model
#'
#' Prognostic feature selection using Cox PH model with regularization methods
#' (Ridge with alpha=0; Lasso with alpha=1;ElasticNet with alpha in the range 0,1)
#'
#' @param feature_data A data.frame contains tumors in rows and aggregated cellular features in columns;
#' row.names represent the tumor identifiers
#' @param surv_data A data.frame contains tumors in rows and survival data in columns:
#' "cens" for censoring status; "time" for surviving time; row.names represent the tumor identifiers
#' @param log10_featureNms A vector of characters of feature names for log10 tranformation prior to z-score
#' transformation; all other features are transformed to z-scores
#' @param output_dir A file.path indicating the output directory
#' @param no_of_runs An integer of the number of repeats of surviving test
#' @param regularization_alpha A fraction (0 to 1) representing the alpha value
#' for \link[glmnet]{cv.glmnet}
#' @return A data.frame contains tumor-level feature summaries saved in
#'   \code{output_dir}
#' @importFrom survival Surv
#' @import glmnet
#' @importFrom stats model.matrix coef as.formula
#' @export
featureSelection_ProgSigFeature <- function (feature_data=NULL, surv_data=NULL, log10_featureNms=NULL,
                                       regularization_alpha=1,
                                          output_dir = NULL, no_of_runs=10){
  ## ======================
  ## exclude features with many 0's or NA's
  ## ======================
  colids_rm <- c()
  fnm <- file.path(output_dir, "FS_PrognosisSignificance_log.txt")
  fnm_con <- file(fnm, open="wt")

  sink(fnm_con)
  sink(fnm_con, type="message")
  dim(feature_data)
  cat('\n\n')

  for(ii in 1:ncol(feature_data)){
    zz <- sum(feature_data[,ii]==0, na.rm = TRUE)
    zzz <- sum(is.na(feature_data[,ii]==0))
    zzzz <- sd(feature_data[,ii] , na.rm = TRUE)

    if(zz > (nrow(feature_data) * 0.80)){
      cat(colnames(feature_data)[ii],': no. of 0s=',zz,'\n')
      colids_rm <- c(colids_rm, ii)}
    if(zzz > (nrow(feature_data) * 0.60)){
      cat(colnames(feature_data)[ii],': no. of NAs=',zzz,'\n')
      colids_rm <- c(colids_rm, ii)}
    if(zzzz == 0){
      cat(colnames(feature_data)[ii],': is invariant, stdev=0 \n')
      colids_rm <- c(colids_rm, ii)}
  }
  cat('\n\n')
  if(length(colids_rm)>0) feature_data <- feature_data[, -colids_rm]
  dim(feature_data)

  sink(type="message")
  sink()


  ## ======================
  ## merge surv data and feature data
  ## ======================
  data <- merge(feature_data, surv_data, by='row.names')
  rownames(data) <- data$Row.names
  data$Row.names <- NULL

  ## ======================
  ## remove cases with survival time ==0
  ## ======================
  d <- data[data$time > 0,]

  ## ======================
  ## exclude cases with NA values
  ## ======================
  rs_NA <- rowSums(is.na(d))
  d <- d[rs_NA==0, ]

  ## ======================
  ## data transformation
  ## ======================
  feature_nms <- setdiff(colnames(d),c("cens","time"))
  for (cc in feature_nms){

    if(cc %in% log10_featureNms){
      cat('LOG10:', cc, '\n')
      d[,cc] <-log10_transform(d[,cc], method = 'log10')
    }
    d[,cc] <- scale(d[,cc])
  }

  ## ======================
  ## Cox PH model with lasso selection
  ## ======================
  ns_terms <- paste0("poly(",feature_nms,", 3, raw = TRUE)")
  x <- model.matrix(as.formula(paste( "~",paste0(ns_terms, collapse = "+"))), d)
  y <- Surv(d$time, d$cens)
  all_runs <- NULL

  set.seed(999)
  for (ii in 1:no_of_runs){
    cat('run# ',ii, '\n')

    cv.fit <- cv.glmnet(x, y, family="cox", alpha=regularization_alpha, nfolds=5)

    ## identify the best model
    min_lambda_id <- which(cv.fit$lambda==cv.fit$lambda.min)
    cv.fit$glmnet.fit$beta[1:nrow(cv.fit$glmnet.fit$beta), min_lambda_id]

    ## identify the selected features
    retained_features <- as.data.frame(as.matrix(coef(cv.fit, s = "lambda.min")))
    retained_features$terms <- rownames(retained_features)
    colnames(retained_features)[1] <- paste0('run', ii)

    ## saving intermediate results
    if(is.null(all_runs)) {
      all_runs <- retained_features
    }else{
      all_runs <- merge(all_runs, retained_features, by='terms')
    }

    #save(all_runs, file = file.path(output_dir,paste0('Lasso_temp.Rdata')))
  }# end no_of_runs

  save(all_runs, file = file.path(output_dir,paste0('FS_PrognosisSignificance_',no_of_runs,'runs.Rdata')))
}

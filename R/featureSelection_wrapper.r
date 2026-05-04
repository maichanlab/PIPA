#' Feature selection wrapper function
#'
#' Feature selection wrapper function: step (1) compute_aggregated_feature,
#' step (2) featureSelection_ProgSigFeature, step (3) featureSelection_patientGrouping
#' step (4) featureSelection_extremeProgGps,
#' step (5) make visualizatin plots for selecting the optimal subset of cell features
#' @param data A data.frame contains cells in rows and features in columns;
#' the column 'tumor_ids' represents the tumor identifiers
#' input argument for \code{compute_aggregated_feature}
#' @param log10_featureNms A vector of characters of feature names for log10 tranformation prior to z-score
#' transformation; all other features are transformed to z-scores;
#' input argument for \code{featureSelection_ProgSigFeature}
#' @param no_of_runs An integer of the number of repeats of surviving test;
#' input argument for \code{featureSelection_ProgSigFeature}
#' @param min_cluster_size An integer indicating the minimum cluster size i.e. cluster/tumor subtypes containing
#' fewer tumors than this number will be excluded in the survival analysis
#' input argument for \code{featureSelection_extremeProgGps}
#' @param surv_data A data.frame: tumors in rows (tumor_ids as row.names), columns must include: "time", "cens";
#' "cens" for censoring status; "time" for surviving time; row.names represent the tumor identifiers;
#' if there are additional columns, they are treated as confounders and will not be used for feature selection
#' input argument for \code{featureSelection_ProgSigFeature} & \code{featureSelection_extremeProgGps}
#' @param output_dir A file.path string of output directory
#' @import ggplot2 dplyr
#' @importFrom stats pchisq
#' @export

featureSelection_wrapper <- function(data=NULL, surv_data=NULL,
                                     log10_featureNms=NULL,  no_of_runs=10,
                                     output_dir=NULL,
                                     min_cluster_size= 20){
  avg_concordance = detection_rate = key = value = prog_group = NULL

  ## =========================
  ## (STEP 1) aggregate cell features to tumpr/patient-level
  ## =========================
  compute_aggregated_feature(data=data, output_dir = output_dir)

  ## =========================
  ## (STEP 2) evaluate prognosis signficance of individual aggregate cell features
  ## {an output coefficient >0 indicates significant association}
  ## =========================
  ## prepare surv data
  surv_data <- surv_data[, c("cens","time")]
  ## load aggregate patient feature data
  aggregated_feature_data <- read.csv(file.path(output_dir,'aggregate_patient_features.csv'), as.is = TRUE, row.names = 1)
  ## evaluate prognostic significance
  log10_featureNms <- grep(x=colnames(aggregated_feature_data), pattern = paste0(log10_featureNms, collapse = "|"),
                           value = TRUE)

  featureSelection_ProgSigFeature(feature_data=aggregated_feature_data, surv_data=surv_data,
                                         log10_featureNms=log10_featureNms,
                                         output_dir = output_dir, no_of_runs=no_of_runs)

  ## =========================
  ## (STEP 3) identify patient groups shared similar cell feature characteristics.
  ## Using cell features obtained at each 10% increment of detection rates
  ## (determined by  function "featureSelection_ProgSigFeature"), tumors are divided into
  ## groups which are evaluated for prognostic significance.
  ## =========================
  ## load aggregate patient feature data determined by "featureSelection_ProgSigFeature" function
  selected_feature <- get(load(file = file.path(output_dir,
                                                paste0('FS_PrognosisSignificance_',no_of_runs,'runs.Rdata'))))
  ## find patient groups and evaluate their prognostic associations
  featureSelection_patientGrouping(cell_feature_data=aggregated_feature_data,
                                   surv_data= surv_data,
                                   min_cluster_size = min_cluster_size,
                                   log10_featureNms=log10_featureNms,
                                   plot_bnm='', output_dir = output_dir,
                                   feature_data = selected_feature)

  ## =========================
  ## (STEP 4) at each unique subset of features identified above
  ## (a) identify patient groups with extreme prognoses (for determining the training tumors),
  ## and (b) compute concordance (as a measure of separation) of the identified
  ## patient groups in (a) with all other patient groups
  ## =========================
  ## results container for all candidate feature subsets
  all_res <- NULL
  ## folder names for each candidate feature subsets
  feature_subset_fnms <- list.files(path = output_dir, pattern = "^detection_rate")
  for (ff in feature_subset_fnms){
    output_dir2 <- file.path(output_dir, ff)

    # load patient groups
    tumor_clusters<-read.csv(file = file.path(output_dir2,'FS_patientGrouping.csv'), as.is = TRUE)

    if(length(unique(tumor_clusters$cluster_no))<=1)next
    # identify patient groups with extreme prognoses and compute the corresponding concordance for feature selection
    res <-featureSelection_extremeProgGps(surv_data=surv_data, subtype=tumor_clusters,
                                          output_dir=output_dir2, min_cluster_size= min_cluster_size)
    cat(res$most_separated_curves$better_surv, ' vs ', res$most_separated_curves$worse_surv,
        ': pval=', res$most_separated_curves$pval, ", avg concor=", res$most_separated_curves_meanConcordance,'\n')

    # prepare data for visualization plot
    # load selected features
    selected_features <- read.csv(file = file.path(output_dir2,'selectedAggFeatures.csv'), as.is = TRUE)
    selected_features <- gsub(x=selected_features[,1], pattern = '^stdev.|Q1.|Q2.|Q3.|min.|max.', replacement = '')
    selected_features <- unique(selected_features)
    # combine features & concordance data
    tmp <- data.frame(detection_rate=ff, avg_concordance=res$most_separated_curves_meanConcordance,
                      best_surv_gp=res$most_separated_curves$better_surv,
                      worst_surv_gp=res$most_separated_curves$worse_surv,
                      best_surv_gp_size=res$most_separated_curves$better_surv_size,
                      worst_surv_gp_size=res$most_separated_curves$worse_surv_size,
                      no_of_features=length(selected_features))
    tmp$detection_rate <- gsub(x=tmp$detection_rate, pattern = "detection_rate", replacement = "")
    # consolidate results across candidte feature subsets
    all_res <- rbind(all_res, tmp)
  }
  # exclude cell feature subset contain single feature
  all_res<- all_res[all_res$no_of_features >= 1,]

  ## =========================
  ## (STEP 5) Visualizatin plots for selecting the optimal subset of cell features
  ## =========================
  ## ----------------
  ## barplot of concordance vs. %detection_rate
  ## ----------------
  plot_concordance <- ggplot(all_res, aes(y=avg_concordance, x=detection_rate)) +
    geom_bar(stat="identity",position=position_dodge())+
    theme_bw()+
    theme(axis.text = element_text(color='black',size = 7),
          axis.title = element_text(color='black',size=7),
          panel.border = element_rect(color='black',size=0.5),
          panel.grid = element_blank(),
          legend.spacing.x = unit(0.2, 'cm'),
          legend.title =  element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.position = 'bottom',
          legend.direction = 'horizontal')+
    ylab("Mean concordance")+ xlab('Dection rate')+
    guides(fill = guide_legend(title = 'phenotype', title.position = "left"))+
    guides(fill = guide_legend(override.aes = list(size = 3)))

  plot_fnm <- file.path(output_dir, "plot_concordance.pdf")
  ggsave(filename=plot_fnm, plot=plot_concordance)


  ## ----------------
  ## barplot of cluster size vs. %detection_rate
  ## ----------------
  gp_size <- all_res[,c("detection_rate","best_surv_gp_size", "worst_surv_gp_size")]
  cols <- grep(x=colnames(gp_size), pattern = "gp_size")
  gp_size <- gp_size%>% tidyr::gather(key,value, all_of(cols))
  gp_size$prog_group <- ifelse(grepl(x=gp_size$key, pattern = "best_surv"), "best surv", "worst surv")
  plot_clusterSize <- ggplot(gp_size, aes(y=value, x=detection_rate, fill=prog_group)) +
    geom_bar(stat="identity",position=position_dodge())+
    theme_bw()+
    theme(axis.text = element_text(color='black',size = 7),
          axis.title = element_text(color='black',size=7),
          panel.border = element_rect(color='black',size=0.5),
          panel.grid = element_blank(),
          legend.spacing.x = unit(0.2, 'cm'),
          legend.title =  element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.position = 'right',
          legend.direction = 'vertical')+
    ylab("Cluster size")+ xlab('Detection rate')+
    guides(fill = guide_legend(title = 'Prognostic group', title.position = "top"))


  plot_fnm <- file.path(output_dir, "plot_clusterSize.pdf")
  ggsave(filename=plot_fnm, plot=plot_clusterSize)
}

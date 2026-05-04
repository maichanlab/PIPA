#' Feature selection - identify and evaluate patient groups with extreme prognoses
#'
#' Feature selection - identify and evaluate patient groups with extreme prognoses.
#' (1) identifying patient groups (among all output from \code{featureSelection_patientGrouping}) with extreme prognoses
#' i.e., most unfavorable vs. favorable
#' (2) computing concordance of Kaplan-Meier estimates as a measure of noise & separation
#'
#' ## step (2) at each unique subset of features tested in step(1),
## (a) identify patient groups with extreme prognoses (i.e. training tumors),
## and (b) compute concordance (as a measure of separation) of the identified
## patient groups with all other patient groups

#' @param min_cluster_size An integer indicating the minimum cluster size i.e. cluster/tumor subtypes containing
#' fewer tumors than this number will be excluded in the survival analysis
#' @param surv_data A data.frame: tumors in rows (tumor_ids in row.names), columns must include: "time", "cens";
#' if there are additional columns, they are treated as confounders
#' @param out_bnm A character string for naming the output result in .RData
#' @param subtype A data.frame representing the tumor subtype data; it contains 2 columns named as tumor_ids
#' and cluster_no
#' @param save_res A boolean indicating if the results to be saved
#' @param output_dir A file.path string of output directory
#' @importFrom stats pchisq
#' @export

featureSelection_extremeProgGps <- function(surv_data=NULL, subtype=NULL,
                                                  save_res=TRUE,
                                                  output_dir=NULL,
                                                  out_bnm = 'FS_extremeProgGps',
                                                  min_cluster_size= 20){

  ## ======================
  ## complete cases with survival
  ## ======================
  surv_data <- surv_data[complete.cases(surv_data[,c("cens","time")]),]


  #### merge surv_data & cluster results
  m <- merge(subtype, surv_data, by.x='tumor_ids', by.y = 'row.names')
  cluster_no <- NULL
  m$cluster_no <- factor(m$cluster_no)
  no_cluster <- length(unique(m$cluster_no))


  #### exclude cluster with small size
  freq <- data.frame(table(m$cluster_no))
  if(max(freq$Freq) < min_cluster_size){
    warning('ALL clusters are smaller than the minimum required size!\n')
    return(NULL)
  }
  outlier_clusterids <- freq$Var1[which(freq$Freq < min_cluster_size)]
  m <- m[!m$cluster_no %in% outlier_clusterids,]
  m$cluster_no <- droplevels(m$cluster_no)
  if( length(unique(m$cluster_no))<=1 ){
    warning('Only 1 cluster is found!\n')
    return(NULL)
  }
  rm(freq)


  #### factorize cluster id
  m$cluster_no <- as.factor(m$cluster_no)
  all_res <- NULL
  best_res <- data.frame(better_surv=-1, worse_surv=-1, concordance=0, pval=100)
  cluster_ids <- levels(m$cluster_no)
  all_concordance <- 0
  for(c1 in 1:length(cluster_ids)){
    for(c2 in (c1+1):length(cluster_ids)){
      if(c2 > length(cluster_ids))break
      m_sub <- m[m$cluster_no%in% c(cluster_ids[c1],cluster_ids[c2]),]
      m_sub$cluster_no <- droplevels(m_sub$cluster_no)
      sc <- survConcordance(Surv(time, cens) ~cluster_no, data=m_sub)
      #cat(cluster_ids[c1], " vs ", cluster_ids[c2], '\n:')
      #cat(sc$concordance, '\t')
      all_concordance <- c(all_concordance, sc$concordance)
      surv_diff <- survdiff(Surv(time, cens) ~ cluster_no, data = m_sub)
      p_val <- pchisq(surv_diff$chisq, length(surv_diff$n)-1, lower.tail = FALSE)


      better_surv_gp <-names(surv_diff$n)[surv_diff$obs/surv_diff$exp==min(surv_diff$obs/surv_diff$exp)]
      better_surv_gp <-gsub(x=better_surv_gp, pattern = "cluster_no=", replacement = "")

      worse_surv_gp <- setdiff(c(cluster_ids[c1], cluster_ids[c2]), better_surv_gp)
      tmp <- data.frame(better_surv=better_surv_gp, worse_surv=worse_surv_gp,
                        concordance=sc$concordance, pval=p_val)
      all_res <- rbind(all_res, tmp)

      #cat(p_val, '\n')
      if(best_res$pval > p_val){
        best_res$pval = p_val
        #best_res$c1=cluster_ids[c1]
        #best_res$c2=cluster_ids[c2]
        best_res$better_surv <- better_surv_gp
        best_res$worse_surv <- worse_surv_gp
        best_res$concordance=sc$concordance
      }
    }
  }

  all_res$concordance
  r1 <- which(all_res$better_surv %in% c(best_res$better_surv, best_res$worse_surv))
  r2 <- which(all_res$worse_surv %in% c(best_res$better_surv, best_res$worse_surv))
  most_separated_curves_meanConcordance = mean(all_res$concordance[c(r1,r2)])

  cluster_size <- m %>% group_by(cluster_no) %>% count
  best_res$better_surv_size <- as.integer(cluster_size[cluster_size$cluster_no == best_res$better_surv, "n"])
  best_res$worse_surv_size <- as.integer(cluster_size[cluster_size$cluster_no == best_res$worse_surv, "n"])

  if(save_res)save(all_ccomparisons=all_res, most_separated_curves=best_res,
                   file=file.path(output_dir, paste0(out_bnm, '.RData')))
  return(list(most_separated_curves_meanConcordance=most_separated_curves_meanConcordance,
              most_separated_curves=best_res))
}

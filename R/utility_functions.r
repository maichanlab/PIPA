#' Data transformation function
#'
#' Data transformation function: log10; log(value <= 0) assigned to log(min. value >0)
#' @param data Input data: samples (cells or tumors) in rows, features in columns
#' @param method A character string indicate the type of data transformation
#' (current version does only 1og10 transformation)
log10_transform <- function(data=NULL, method='log10'){
  if(method=='log10'){
    no_NA <- sum(is.na(data))
    if(no_NA>0)warning(no_NA, ' NAs is detected\n')

    z1<-log(data, base = 10)
    z1[is.infinite(z1)] <- min(z1[!is.infinite(z1)])
    return(z1)
  }

}# end log10_transform


#' Survival analysis wrapper function
#'
#' Survival analysis wrapper function
#' @param min_cluster_size An integer indicating the minimum cluster size i.e. cluster/tumor subtypes containing
#' fewer tumors than this number will be excluded in the survival analysis
#' @param KM_xlab A character string for labelling the x-axis in Kaplan-Meier curves
#' @param surv_data A data.frame: tumors in rows (tumor_ids in row.names), columns must include: "time", "cens";
#' if there are additional columns, they are treated as confounders
#' @param plot_width An integer of the PDF file width
#' @param out_bnm A character string for naming the output PDF file
#' @param subtype A data.frame representing the tumor subtype data; it contains 2 columns named as tumor_ids
#' and cluster_no
#' @param user_ref A character indicating the reference group (ie. tumor subtype)
#' @param output_dir A file.path string of output directory
surv_analysis_fn <- function(surv_data=NULL, subtype=NULL, user_ref=NULL,
                             output_dir=NULL,plot_width=8,
                             out_bnm = 'CoxPH',
                             min_cluster_size= 20, KM_xlab='Survivals'){

  ## ======================
  ## identify confounders
  ## ======================
  confounders <- setdiff(colnames(surv_data), c('cens','time'))
  ## ---------------
  ## color palette
  ## ---------------
  col_vec <- c(brewer.pal(9, 'Set1')[-6], brewer.pal(8, 'Set2'), brewer.pal(12,'Set3')[-c(2,12)],
               brewer.pal(8, 'Accent')[-4],brewer.pal(12, 'Paired')[-11])
  ## ======================
  ## complete cases with survival
  ## ======================
  surv_data <- surv_data[complete.cases(surv_data[,c("cens","time")]),]


  #### merge surv_data & cluster results
  m <- merge(subtype, surv_data, by.x='tumor_ids', by.y = 'row.names')
  cluster_no <- NULL
  m$cluster_no <- factor(m$cluster_no)
  no_cluster <- length(unique(m$cluster_no))
  col_vec2 <- col_vec[c(1:no_cluster)]
  names(col_vec2) <- as.character(levels(m$cluster_no))

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


  #### factorize TIPC cluster id
  m$cluster_no <- as.factor(m$cluster_no)
  ## output data container var
  output_HR_CI <- list(); ll <- 1
  #### set individual clusters as reference cluster
  if(is.null(user_ref)){
    for (ref_cluster_no in unique(m$cluster_no)){

      m$cluster_no  <- relevel(m$cluster_no ,ref=as.character(ref_cluster_no))


      ## univariate Cox PH
      if(length(confounders)>0){
        formula_mult <- as.formula(paste("Surv(time, cens)~cluster_no+",paste0(confounders,collapse = '+')))
        res.cox <- survival::coxph(formula_mult, data =  m)
        fileNm <- file.path(output_dir, paste0('multiCox_ref',ref_cluster_no,'.txt'))
        capture.output(summary(res.cox), file = fileNm)
      }

      formula <- as.formula(paste("Surv(time, cens)~cluster_no"))
      res.cox <- survival::coxph(formula, data =  m)
      fileNm <- file.path(output_dir, paste0('uniCox_ref',ref_cluster_no,'.txt'))
      capture.output(summary(res.cox), file = fileNm)
      ## output saving
      output_HR_CI[[ll]] <- summary(res.cox)$conf.int
      names(output_HR_CI)[ll] <- ref_cluster_no
      ll <- ll + 1
    }#end ref_cluster_no
  }else{
    ref_cluster_no <- user_ref
    m$cluster_no  <- relevel(m$cluster_no ,ref=as.character(ref_cluster_no))
    ## univariate Cox PH
    if(length(confounders)>0){
      formula_mult <- as.formula(paste("Surv(time, cens)~cluster_no+",paste0(confounders,collapse = '+')))
      res.cox <- survival::coxph(formula_mult, data =  m)
      fileNm <- file.path(output_dir, paste0('multiCox_ref',ref_cluster_no,'.txt'))
      capture.output(summary(res.cox), file = fileNm)
    }

    formula <- as.formula(paste("Surv(time, cens)~cluster_no"))
    res.cox <- survival::coxph(formula, data =  m)
    fileNm <- file.path(output_dir, paste0('uniCox_ref',ref_cluster_no,'.txt'))
    capture.output(summary(res.cox), file = fileNm)

    ## output saving
    output_HR_CI[[ll]] <- summary(res.cox)$conf.int
    names(output_HR_CI)[ll] <- ref_cluster_no
    ll <- ll + 1
  }#end ref

  ## Kaplan meier curves
  cluster_size <- m %>% group_by(cluster_no) %>% count()
  if(is.null(user_ref))ref_cluster_no <- cluster_size$cluster_no[which.max(cluster_size$n)]

  m$cluster_no  <- relevel(m$cluster_no ,ref=as.character(ref_cluster_no))
  #m$cluster_no <- relevel(m$cluster_no, ref="NST")
  col_vec2 <- col_vec2[match(as.character(levels(m$cluster_no)),names(col_vec2))]

  formula <- as.formula(paste("Surv(time, cens)~cluster_no"))
  fit <- surv_fit(formula, data = m)

  table_ht <- min(length(unique(m$cluster_no))*0.025, 0.7)
  table_ht <- max(table_ht, 0.3)
  pl<-ggsurvplot(fit,ylim = c(0, 1),
                 palette = as.vector(col_vec2),
                 pval = TRUE,
                 pval.size=8,
                 fontsize = c(6),
                 font.title = c(16, "bold"),
                 font.subtitle = c(15, "bold"),
                 font.caption = c(16, "plain"),
                 font.x = c(18, "bold"),
                 font.y = c(18, "bold"),
                 font.tickslab = c(18),
                 font.legend = c(16),
                 risk.table.y.text = FALSE,
                 #trim.strata.names = TRUE,
                 xlab=KM_xlab,
                 risk.table = TRUE, # Add risk table
                 risk.table.col = "strata", # Change risk table color by groups
                 #linetype = "strata", # Change line type by groups
                 tables.height = table_ht,
                 conf.int = FALSE,
                 surv.median.line = "hv", # Specify median survival
                 ggtheme = theme_bw()) # Change ggplot2 theme

  ht <- max(table_ht * 30, 8)
  #ht <- max(length(unique(m$cluster_no))*0.8, 8)

  pl$plot <- pl$plot +
    theme(legend.text = element_text(size = 14, color = "black", face = "bold"))+
    guides(color = guide_legend(size=10,ncol=3,byrow=TRUE))

  plot_fileNm <- file.path(output_dir, paste0('KM_ref',ref_cluster_no,'.pdf'))
  pdf(plot_fileNm, onefile = FALSE, height = ht, width = plot_width)
  print(pl)
  dev.off()

  save(output_HR_CI, file=file.path(output_dir, paste0(out_bnm, '_HR_CI.RData')))
}
#' Consensus matrix updating function
#'
#' Consensus matrix updating function
#' @param clusterAssignments A named vector of length X-samples containing the input cluster assignment of each sample
#' @param m A data.frame of size X-by-X containing the cumulated consensus cluster assignment
#' @param sampleKey A vector of length X containing the sample identifiers
#' @export
connectivityMatrix <- function (clusterAssignments, m, sampleKey) {
  names(clusterAssignments) <- sampleKey
  cls <- lapply(unique(clusterAssignments), function(i) as.numeric(names(clusterAssignments[clusterAssignments %in%
                                                                                              i])))
  for (i in 1:length(cls)) {
    nelts <- 1:ncol(m)
    cl <- as.numeric(nelts %in% cls[[i]])
    updt <- outer(cl, cl)
    m <- m + updt
  }
  return(m)
}# end connectivityMatrix

#' PhenoGraph wrapper function
#'
#' PhenoGraph wrapper function
#' @param data A data frame: samples in rows, features in columns
#' @param seed An integer of seed for reproducibility
#' @param k An integer of the number of neighbors
#' @return A data frame which appends the cluster assignment to input data \code{data}
#' @importFrom stats as.dist complete.cases cutree hclust
#' @import Rphenograph
#' @export

phenograph_clustering <- function(data=NULL, k=30, seed=999){
  data <- data[complete.cases(data), ,drop=FALSE]
  total_cell <- nrow(data)
  set.seed(seed)
  data <- data[sample(x=c(1:total_cell), size = total_cell),,drop=FALSE]

  Rphenograph_out <- Rphenograph(data, k = k)
  missing_case <- setdiff(1:nrow(data),names(igraph::membership(Rphenograph_out[[2]])))
  if(length(missing_case) >0)data <- data[-missing_case, ]
  data$cluster <- factor(igraph::membership(Rphenograph_out[[2]]))

  return(data)}# end phenograph_clustering

#' Cell consensus clustering function
#'
#' Cell consensus clustering function using PhenoGraph to group cells
#' @param data Input data: samples in rows, features in columns; cell_ids in row.names
#' @param batch_size An integer indicating the cell window size i.e. number of cells in each batch run
#' @param consensus_rep An integer indicating the number of repeats for consensus clustering (per window/batch run)
#' @param max_k An integer indicating the maximum number of clusters used in hierachical clustering (i.e. k in cutree)
#' of the consensus matrix, at each bath run
#' @param output_dir A file.path string of output directory
#' @param phenograph_k An integer indicating the number of neighbors in PhenoGraph
#' @param max_batch An integer indicating the maximum number of windows/batchs
#' @importFrom stats as.dist complete.cases cutree hclust

CCC <- function(data=NULL,batch_size=3000,consensus_rep=5, max_k=20,
                output_dir=getwd(), phenograph_k=100,
                max_batch=NULL){

  if(nrow(data) < batch_size) batch_size <- nrow(data)
  ## ======================
  ## randomize global cell ordering
  ## ======================

  set.seed(999)
  rand_order <- sample(x = 1:nrow(data), size = nrow(data))
  rand_data <- data[rand_order, ,drop=FALSE]
  if(is.null(max_batch))  max_batch <- ceiling(nrow(rand_data)/batch_size) #nrow(rand_data) %/% batch_size
  ## generate seeds for window selection
  window_seeds <- sample.int(n=10000, size=max_batch)
  ## generate seeds for PhenoGraph clustering
  PhenoGraph_seeds <- sample.int(n=10000, size=consensus_rep)
  ## ======================
  ## ======================
  ## processing randomly defined cell windows with the same size
  ## ======================
  ## ======================

  batch_count <- 1
  cluster_res <- list()
  ll <- 1
  while(batch_count <= max_batch){

    cat('\n@batch#=',batch_count, ' ( total batches=',max_batch,')\n')

    ## random selection of cells
    set.seed(window_seeds[batch_count])
    data_subset <- rand_data[sample(x = 1:nrow(rand_data), size = batch_size), ,drop=FALSE]

    ############################
    ## consensus clustering preparation: empty matrix
    ############################
    mConsist <- matrix(data = 0, nrow = nrow(data_subset), ncol= nrow(data_subset))
    rownames(mConsist) <- rownames(data_subset)
    colnames(mConsist) <- rownames(data_subset)

    ############################
    ## repeat clustering for consensus_rep
    ############################
    #set.seed(999)

    for (rep in 1:consensus_rep){
      cat('\n\nrepeat#',rep, '\n')
      ## calling PhenoGraph
      pheno_clusters_rep <- phenograph_clustering(data = data_subset, k=phenograph_k, seed = PhenoGraph_seeds[rep])

      ## clustering assignment
      this_assignment <- as.character(pheno_clusters_rep$cluster)
      names(this_assignment) <- rownames(pheno_clusters_rep)
      ## match ordering in mConsist
      this_assignment <- this_assignment[match(rownames(mConsist), names(this_assignment))]
      stopifnot(identical(rownames(mConsist), names(this_assignment)))


      ## write to consensus matrix
      if(rep==1)consensus_mat <- mConsist
      consensus_mat <- connectivityMatrix(this_assignment, consensus_mat,
                                          c(1:ncol(mConsist)))
    }# end consensus_rep


    ## ======================
    ## hierarchical clustering of consensus clusters
    ## ======================
    finalLinkage = "average"
    hc = hclust(as.dist(consensus_rep - consensus_mat), method = finalLinkage)

    ## cutree at max_k
    ct = cutree(hc, k=max_k)

    ## ======================
    ## identify final cluster assignment
    ## ======================
    pheno_clusters <- data.frame(cell_ids=names(ct), phenotype=ct)
    #head(pheno_clusters)

    cluster_res[[ll]] <- pheno_clusters
    names(cluster_res)[ll] <- batch_count
    ll <- ll+1

    ## ======================
    ## saving temporary results
    ## ======================
    #save(cluster_res, consensus_mat, file=file.path(output_dir, paste0('batch',batch_count,'_CCC.RData')))
    #save(cluster_res, consensus_mat, file=file.path(output_dir, paste0('batchX_CCC.RData')))
    batch_count <- batch_count+1
  }# end while

  ############################
  ## save final results
  ############################
  save(cluster_res, file=file.path(output_dir, paste0('cellCluster_byBatches.RData')))
}#end function CCC



#' Meta-clustering of the cell-consensus-clusters
#'
#' Meta-clustering of the cell-consensus-clusters obtained from overlapping data windows
#' @param data A list whose elments are clustering assignment data.frame output from \code{CCC}
merge_CCC <- function(data=NULL){
  ## merge every contingous matrices in the list
  meta_ccc <- list()
  ll <- 1
  ids <- c(1:length(data))
  for (ii in ids){
    if(ii %% 2  == 0) next
    if((ii+1) > length(data)) {
      ## get data
      d1 <- data[[ii]]
      d2 <- meta_ccc[[1]]

    }else{
      ## get data
      d1 <- data[[ii]]
      d2 <- data[[ii+1]]
    }
    cat('\nwindow#',ii, '\n')

    ## merging
    resXY <- merge(d1, d2, by='cell_ids', suffixes = c('.resX', '.resY'))
    if(nrow(resXY)==0) next
    rownames(resXY) <- resXY$cell_ids
    resXY$cell_ids <- NULL
    #head(resXY)
    freq <- as.matrix(table(resXY))
    ## checking overlapping cases
    print(freq)

    ## matching the most enriched cluster: in X
    row_max <- apply(freq, MARGIN = 1, which.max)
    col_max <- apply(freq, MARGIN = 2, which.max)

    ## check if row-wise & column-wise agree
    row_max <- data.frame(rowClust=names(row_max), colClust=colnames(freq)[row_max])
    col_max <- data.frame(rowClust=rownames(freq)[col_max], colClust=names(col_max))

    # col+row cluster ids
    row_max$rc_clust <- paste0(row_max$rowClust,'_', row_max$colClust)
    col_max$rc_clust <- paste0(col_max$rowClust,'_', col_max$colClust)

    ## checking overlapping cases
    consistent_clusterids <- intersect(row_max$rc_clust, col_max$rc_clust)
    consistent_clusterXids <- sapply(strsplit(x = consistent_clusterids, split = '_'), "[[", 1)
    consistent_clusterYids <- sapply(strsplit(x = consistent_clusterids, split = '_'), "[[", 2)
    consistent_clusterXids <- as.integer(consistent_clusterXids)
    consistent_clusterYids <- as.integer(consistent_clusterYids)

    ## clean up cluster in d1 & d2
    d1 <- d1[d1$phenotype%in% (consistent_clusterXids), ]
    d2 <- d2[d2$phenotype%in% (consistent_clusterYids), ]
    # rename phenotype ids in d2
    table(d2$phenotype)
    d2$new_cluster <- 0
    for (YY in c(1:length(consistent_clusterYids))){
      cat(consistent_clusterYids[YY], ' replaced by ', consistent_clusterXids[YY], '\n')

      d2$new_cluster[d2$phenotype == consistent_clusterYids[YY]] <- consistent_clusterXids[YY]
    }
    table(d2$new_cluster)
    d2$phenotype <- d2$new_cluster
    d2$new_cluster <- NULL

    ## merge d1 and d2
    meta_cluster <- rbind(d1,d2)
    meta_cluster <- meta_cluster[!duplicated(meta_cluster$cell_ids), ]
    sum(duplicated(meta_cluster$cell_ids))
    table(meta_cluster$phenotype)

    meta_ccc[[ll]] <- meta_cluster
    ll <- ll +1
  }# end pairwise merging
  return(meta_ccc)
}# end function merge_CCC

#' Forest plot summarizing Cox PH regression analysis results
#'
#' Forest plot summarizing Cox PH regression analysis results - supporting function for \code{ForesPlot_denQ_CoxPH}
#' @param HR_CI A data.frame containing the HR and CI data
#' @param param_bnm A character string of the identifier of the experimental setting of interest
#' to be used in \link[ggplot2]{facet_wrap}; and it will appended to the PDF output filename
#' @param tissue_areas A vector of characters indicating the tissue regions in which density quartiles
#' are computed are to be used for making forest plots; it has to be consistent with the tissue region
#' naming found under the 'area' column in input data HR_CI
#' @param output_dir A file.path string of output directory
#' @return Number PDF files equal to the length of tissue_areas
#' @param plot_width An integer to be passed to the 'width' argument in
#' \link[ggplot2]{ggsave}
#' @param plot_height An integer to be passed to the 'height' argument in
#' \link[ggplot2]{ggsave}
#' @import dplyr ggplot2
forest_plot   <- function(HR_CI=NULL, output_dir= NULL,
                          param_bnm = 'test',tissue_areas = c('S', 'T'),
                          plot_width=8, plot_height=8){
  HR<-Lasso_selection<- cat_model <- head<- lower_CI <-model<- upper_CI<- NULL
  ## ===============
  ## ordering
  ## ===============
  HR_CI$cat <- factor(HR_CI$cat, levels = rev(c("Q1","Q2","Q3","Q4")))
  HR_CI$model <- factor(HR_CI$model, levels = (c("univariable","multivariable")))
  HR_CI$phenotype <- factor(HR_CI$phenotype, levels = (c("PP","IP","GP")))

  ## ===============
  ## convert to wide format: separate columns for HR and CI
  ## ===============
  HR_CI$cat_model <- paste0(HR_CI$cat, '_', HR_CI$model)
  head(HR_CI)

  HR_CI$cat_model <- factor(HR_CI$cat_model,
                            levels = rev(c("Q1_univariable","Q1_multivariable",
                                           "Q2_univariable","Q2_multivariable",
                                           "Q3_univariable","Q3_multivariable",
                                           "Q4_univariable","Q4_multivariable")))

  HR_CI$cat_model
  HR_CI$param_bnm <- HR_CI[,param_bnm]

  ## ===============
  ## loop over tissue area
  ## ===============

  for( tt in tissue_areas){
    tmp <- HR_CI[HR_CI$area == tt, ]

    ## remove non-ref quartiles not selected
    tmp$Lasso_selection <- 'yes'
    tmp$Lasso_selection[tmp$HR==1 & tmp$cat != 'Q1'] <- 'no'
    table(tmp$Lasso_selection)
    size_vec <- c('yes'=0.5, 'no'=0)

    p = ggplot(data=tmp,
               aes(x = cat_model,y = HR, ymin = lower_CI, ymax = upper_CI , shape=cat))+
      geom_pointrange(aes(col=model, size=Lasso_selection))+
      scale_size_manual(name = "Lasso selection",values = size_vec)+
      geom_hline(aes(fill=model),yintercept =1, linetype=2)+
      xlab('')+
      ylab("Hazard Ratio (95% Confidence Interval)")+
      geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI,col=model),width=0.3,cex=0.5)+
      facet_grid(param_bnm ~ phenotype, scales = 'free_x')+
      theme_bw()+
      theme(#plot.title=element_text(size=16),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=7, color='black'),
        axis.title=element_text(size=7),
        panel.border = element_rect(color = 'black', size = 0.3),
        strip.background  = element_rect(color = 'black', size = 0.3, fill="white"),
        panel.grid = element_blank(),
        strip.text = element_text(size = 7),
        #legend.margin = margin(0,0,0,1, unit="cm"),
        legend.position = 'bottom', legend.direction = 'horizontal',
        legend.box = 'horizontal ',
        legend.title=element_text(size=7),legend.text = element_text(size=7))+
      labs(shape='Density quartiles') +
      scale_color_discrete(name = "Cox model",
                           labels = c("univariate",
                                      "multivariate"))+
      coord_flip() +
      guides(color=guide_legend(nrow=1,byrow=TRUE, override.aes = list(size = 0.5)))+
      guides(fill=guide_legend(nrow=1,byrow=TRUE, override.aes = list(size = 0.5)))+
      guides(shape=guide_legend(nrow=1,byrow=TRUE, override.aes = list(size = 0.5)))


    setwd(output_dir)
    ggsave(plot = p, width = plot_width, height = plot_height,
           filename = paste0(tt, '_forest_plot.pdf'))
  }

}#end function forest_plot


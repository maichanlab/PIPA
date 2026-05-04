#' Construct meta_cluster
#'
#' Construct meta_cluster using training cells
#'
#' @param cell_data A data.frame contains cells in rows and features in columns; cell_ids in row.names
#' the column 'tumor_ids' are the tumor identifiers; the column 'core_ids' are the core identifiers;
#' the column 'cell_ids' are the cell identifiers containing core_ids+coordinates
#' @param log10_featureNms A vector of characters of feature names whose values to be log10-transformed prior to scaling (z-score)
#' @param consensus_rep An integer indicating the number of repeats for consensus clustering (per window/batch run);
#' to be passed to \code{CCC}
#' @param prognosis_delineation A boolean indicating if the resultant phenotypes to be delineated by tumor prognosis
#' subtypes (as identified in Step 1)
#' @param cell_batch_clustering A boolean: TRUE (Step 2) cell clustering in batches will be performed; FALSE
#' previous results will be loaded for merging
#' @param cell_batch_merging A boolean: TRUE (Step 2) the merging of cell clusters from individual batches
#' will be performed;
#' @param User_ClusterNo An integer: User-specified number of meta clusters

#' @param phenograph_k An integer indicating the number of neighbors in PhenoGraph used in \code{CCC}
#' @param max_batch An integer indicating the maximum number of windows/batchs used in \code{CCC}
#' @param batch_size An integer indicating the cell window size i.e. number of cells in each batch run;
#' to be passed to \code{CCC}
#' @param no_of_random_merging An integer indicating the number merging with different (randomized) ordering of
#' cell windows; for identifying the inherent no. of phenotypes
#' @param input_dir A file.path string of input directory
#' @param output_dir A file.path string of output directory that the resultant phenotyping model is saved
#' @param seed An integer for reproducibility

#' @param remove_cor_features A boolean indicating if highly correlated features should be removed
#' @export
#' @import ggplot2 RColorBrewer
#' @importFrom stats cor
#' @importFrom corrplot corrplot
#' @importFrom caret findCorrelation
#' @importFrom amap Dist
#'
meta_cluster <- function(cell_data=NULL,log10_featureNms=NULL, batch_size=3000,
                                   input_dir=NULL, output_dir=NULL, prognosis_delineation=TRUE,
                                   cell_batch_clustering=TRUE, cell_batch_merging=TRUE,
                                   User_ClusterNo=NULL,
                                   max_batch=NULL, phenograph_k=100, consensus_rep=5,
                                   remove_cor_features=FALSE,
                                   no_of_random_merging=100, seed=999){

  if(cell_batch_clustering){
    if(!cell_batch_merging)stop('cell_batch_clustering=TRUE, cell_batch_merging cannot be FALSE\n')
  }
  # if(max_batch < 2){
  #   stop('max_batch has to be more than 1\n')
  # }

  ## ======================
  ## parameters documentation
  ## ======================
  setwd(output_dir)
  capture.output(file = 'meta_cluster.txt',append = FALSE,
                 cat('\n\n',format(Sys.time(), "%a %b %d %X %Y"),
                     '\nprognosis_delineation=',prognosis_delineation,
                     '\nUser_ClusterNo=',User_ClusterNo,
                     '\nno_of_random_merging=' , no_of_random_merging))
  if(cell_batch_clustering){
    capture.output(file = 'meta_cluster.txt',append = TRUE,
                   cat('\nbatch_size=',batch_size, '\nbatch_no=',max_batch,
                       '\nphenograph_k=',phenograph_k, '\nconsensus_rep=',consensus_rep,
                       '\nlog10_featureNms=',paste0(log10_featureNms,collapse = '\t')))
  }


  ## ======================
  ## construct cell/core/tumor identifier data.frame:
  ## ======================
  ID <- cell_data[,c("core_ids","tumor_ids")]
  ID <- unique(ID)
  rownames(cell_data) <- cell_data$cell_ids
  cell_data[,c("cell_ids","core_ids","tumor_ids")] <- NULL

  ## ======================
  ## assign patient groups to best, intermediate, worst surviving categories
  ## ======================
  setwd(input_dir)
  PatientGps <- read.csv('FS_patientGrouping.csv',as.is = TRUE)
  loaded_vars <- load('FS_extremeProgGps.RData')
  input <- list(var1 = get(loaded_vars[1]), var2 = get(loaded_vars[2]))
  PatientGps_pval <- input[[1]]
  bestSurv_gp <- PatientGps_pval$better_surv[which.min(PatientGps_pval$pval)]
  worstSurv_gp <- PatientGps_pval$worse_surv[which.min(PatientGps_pval$pval)]

  PatientGps$prog_gps <- 'intermediate'
  PatientGps$prog_gps[PatientGps$cluster_no %in% bestSurv_gp] <- 'best'
  PatientGps$prog_gps[PatientGps$cluster_no %in% worstSurv_gp] <- 'worst'
  ## ======================
  ## filter prognostic subtypes for best & worst
  ## ======================
  PatientGps <- PatientGps[PatientGps$prog_gps %in% c('best','worst'), ,drop=FALSE]

  cat('Tumor subtyping cluster size (user selection), best vs worst: ',
      table(PatientGps$prog_gps[PatientGps$prog_gps %in% c('best','worst')]), '\n')

  ## ======================
  ## load user selected (cut-off) prognostic features
  ## ======================
  setwd(input_dir)
  selected_features <- read.csv(file = file.path('selectedAggFeatures.csv'), row.names = NULL)

  ## ======================
  ## clean up the (aggregate summary) feature names to unique sub-cellular features
  ## ======================
  selected_features <- gsub(x=selected_features[,1], pattern = '^stdev.|Q1.|Q2.|Q3.|min.|max.', replacement = '')
  selected_features <- unique(selected_features)
  cat('selected features:: ', paste(selected_features, collapse = ', '), '\n')

  ## ======================
  ## filter cell data by selected features
  ## ======================
  cell_data <- cell_data[, colnames(cell_data) %in% selected_features, drop=FALSE]
  stopifnot(ncol(cell_data)==length(selected_features))
  ## ======================
  ## exclude cells not in the best and worst prognosis subtypes
  ## ======================
  ID <- ID[ID$tumor_ids %in% PatientGps$tumor_ids,]
  stopifnot(length(unique(ID$tumor_ids)) == nrow(PatientGps))

  cell_data$core_ids <- sapply(strsplit(x=rownames(cell_data), split = 'XY'), "[[", 1)
  cell_data <- cell_data[cell_data$core_ids %in% ID$core_ids, ]
  cell_data$core_ids <- NULL

  ## ======================
  ## data transformation
  ## ======================
  for (cc in 1:ncol(cell_data)){
    feature_cc <- colnames(cell_data)[cc]
    if(feature_cc %in% log10_featureNms){
      cat('LOG10:', feature_cc, '\n')
      cell_data[,cc] <-log10_transform(cell_data[,cc], method = 'log10')
    }
    cell_data[,cc] <- scale(cell_data[,cc])
  }

  ## ======================
  ## remove invariant (i.e. sd=0) features
  ## ======================
  invariant_feature_id <- which(apply(cell_data, MARGIN = 2, sd)==0)
  if(length(invariant_feature_id)>0){
    cell_data <- cell_data[, -invariant_feature_id]
    write.csv(invariant_feature_id,file = file.path(output_dir, 'invariant_features.csv'))
  }

  ## ======================
  ## remove highly correlated features
  ## ======================
  if(length(selected_features)>=2){
    correlationMatrix <- cor(cell_data)

    highlyCorrelated <- findCorrelation(correlationMatrix)
    if(length(highlyCorrelated)>0 & remove_cor_features){
      cat('Exclude highly correlated features: ', colnames(cell_data)[highlyCorrelated], '\n')
      cell_data <- cell_data[, -highlyCorrelated]
    }
  }
  profFeature_retained <- data.frame(features=colnames(cell_data))
  # write.csv(profFeature_retained,file = file.path(output_dir, 'correlatedFeaturesRemoved.csv'),
  #           row.names = FALSE)
  write.csv(profFeature_retained,file = file.path(input_dir, 'selectedCellFeatures.csv'), row.names = FALSE)
  ## ======================
  ## correlation plot
  ## ======================
  if(remove_cor_features){
    correlationMatrix <- cor(cell_data)
    pdf(file.path(output_dir,'selectedFeatures_clean_corr.pdf'))
    corrplot(correlationMatrix)
    dev.off()
  }

  if(cell_batch_clustering){
    ## ======================
    ## consensus cell clustering in randomly selected cell windows (DIVIDE)
    ## ======================
    #CCC(data=cell_data, output_dir=output_dir, max_batch = 2, consensus_rep = 2)

    CCC(data=cell_data, output_dir=output_dir, batch_size=batch_size,
        consensus_rep=consensus_rep, phenograph_k=phenograph_k, max_batch=max_batch)

  }# end cell_batch_clustering

  ## ======================
  ## cell_batch_merging
  ## ======================
  if(cell_batch_merging){
    ## ======================
    ## merge consensus cell clusters (CONQUER)
    ## ======================
    cluster_res <- get(holder <- load(file=file.path(output_dir, paste0('cellCluster_byBatches.RData'))))

    ## ======================
    ## ======================
    ## repeat merging with different ordering of cell windows
    ## ======================
    ## ======================
    size <- ceiling(length(cluster_res)*1)
    set.seed(999)
    seeds <- sample(x = 100000, size = no_of_random_merging)

    #### effect of ordering of merging
    final_metaClust <- list()
    all_res <- list()
    ll <- 1
    for (r in seeds){
      set.seed(r)
      #rand_order <- sample(x = c(1:length(cluster_res)), size = length(cluster_res))
      rand_order <- sample(x = c(1:length(cluster_res)), size = size)
      cluster_res_temp <- cluster_res[rand_order]

      if(length(cluster_res_temp)<=1){
        meta_clusters <- cluster_res_temp[[1]]
      }else{
        meta_clusters_temp <- merge_CCC(data=cluster_res_temp)
        while (length(meta_clusters_temp)>1){
          meta_clusters_temp <- merge_CCC(data=meta_clusters_temp)
        }

        if(length(meta_clusters_temp) == 0){
          stop("NO overlapping cells across batches! Please increase batch size.")
        }
        meta_clusters<-meta_clusters_temp[[1]]
      }


      final_metaClust[[ll]] <- meta_clusters
      names(final_metaClust)[ll] <- r

      no_final_pheno <- length(unique(meta_clusters$phenotype))
      total_cells_remains <- nrow(meta_clusters)
      res <- data.frame(seed=r, no_final_pheno=no_final_pheno, total_cells_remains=total_cells_remains)
      all_res[[ll]] <- res
      ll <- ll+1
    }
    all_res_df <- do.call(rbind, all_res)
    rm(cluster_res)
    #save(final_metaClust, all_res_df, file= file.path(output_dir,'final_metaClust.RData'))

    ## ======================
    ## check phenotype assignment consistency of runs resulted in max. no. of inherent meta-clusters i.e.phenotypes
    ## ======================
    if(is.null(User_ClusterNo)){
      most_freq_no <- all_res_df %>% group_by(no_final_pheno) %>% count()
      mostFreqClusterNo <- most_freq_no$no_final_pheno[which.max(most_freq_no$n)]
    }else{
      mostFreqClusterNo <- User_ClusterNo
    }

    all_res_df2 <- all_res_df[all_res_df$no_final_pheno == mostFreqClusterNo, ]
    maxClust_seed <- all_res_df2$seed

    consistency_count <- 0
    total_count <- 0
    consistent_seeds <- c()
    for(i in 1:length(maxClust_seed) )  {
      mc <- maxClust_seed[i]
      id <- which(names(final_metaClust)==mc)
      mc_x <-final_metaClust[[id]]

      for(i2 in i:length(maxClust_seed)){
        mmc <- maxClust_seed[i2]
        id2 <- which(names(final_metaClust)==mmc)

        total_count <- total_count+1
        cat(mc,' vs ', mmc, '\n')
        mmc_x <-final_metaClust[[id2]]
        mmm <- merge(mc_x,mmc_x, by='cell_ids')
        freq <- table(mmm[,c(2,3)])

        ## matching the most enriched cluster: in X
        row_max <- apply(freq, MARGIN = 1, which.max)
        col_max <- apply(freq, MARGIN = 2, which.max)
        ## check if row-wise & column-wise agree
        row_max <- data.frame(rowClust=names(row_max), colClust=colnames(freq)[row_max])
        col_max <- data.frame(rowClust=rownames(freq)[col_max], colClust=names(col_max))
        # col+row cluster ids
        row_max$rc_clust <- paste0(row_max$rowClust,'_', row_max$colClust)
        col_max$rc_clust <- paste0(col_max$rowClust,'_', col_max$colClust)

        if(sum(row_max$rc_clust %in% col_max$rc_clust) == mostFreqClusterNo){
          cat(mc,' & ', mmc, ' results are consistent\n')
          consistency_count <- consistency_count+1
          consistent_seeds <- c(consistent_seeds, mc, mmc)
        }
      }
    }#end maxClust_seed

    ## ======================
    ## identify the run with highest consistency with other runs
    ## ======================
    #save(consistent_seeds, file= file.path(output_dir,'consistent_seeds.RData'))

    consistent_seeds <- data.frame(consistent_seeds)
    consistent_seeds <- consistent_seeds %>% group_by(consistent_seeds)%>% count()
    opt_seed <- consistent_seeds$consistent_seeds[which.max(consistent_seeds$n)]
    opt_seed <- as.character(opt_seed)
    meta_clusters <- final_metaClust[[opt_seed]]

    ## ======================
    ## saving the optimal training model/ meta clusters
    ## ======================
    cat('Meta cluster size : ',table(meta_clusters$phenotype), '\n')
    save(meta_clusters, file = file.path(output_dir, 'raw_meta_clusters.RData'))

    ## ======================
    ## plot: distribution of the no. of inherent phenotypes
    ## ======================
    all_res_df$no_final_pheno <- as.factor(all_res_df$no_final_pheno)
    pl <- ggplot(data = all_res_df, aes(x=no_final_pheno, total_cells_remains, fill=no_final_pheno))+
      geom_boxplot() + geom_jitter(color='lightgray')+ theme_bw()

    plot_fnm <- paste0('no_of_phenotypes_across',no_of_random_merging,'runs.pdf')
    setwd(output_dir)
    ggsave(plot = pl, filename = plot_fnm)

    ## ======================
    ## no. of phenotypes distribution
    ## ======================
    freq <- as.matrix(table(all_res_df$no_final_pheno))
    freq <- data.frame(no_of_phenotype =rownames(freq), no_of_runs=freq[,1])
    fractions <- prop.table(table(all_res_df$no_final_pheno))
    freq$fraction_of_runs <- fractions
    write.csv(freq, file = file.path(output_dir,'no_of_phenotype_byRuns.csv'), row.names = FALSE)
    ## ======================
    ## get core and tumor ids on meta-clusters
    ## ======================
    meta_clusters$core_ids <- sapply(strsplit(x= as.character(rownames(meta_clusters)), split = 'XY'), "[[", 1)

    #head(meta_clusters)

    ## ======================
    ## delineate meta clusters into good and bad phenotypes
    ## ======================
    if(prognosis_delineation){
      meta_clusters <- merge(meta_clusters, ID[, c('core_ids','tumor_ids')], by = 'core_ids', all.x = TRUE)
      meta_clusters <- merge(meta_clusters, PatientGps[,c("tumor_ids","prog_gps"), drop=FALSE],
                             by = 'tumor_ids')
      ## filter for cells in the best and worst tumor subtypes
      #meta_clusters <- meta_clusters[meta_clusters$prog_gps %in% c('best','worst'),]

      ## ======================
      ## checking
      ## ======================
      meta_clusters_goodCells <- meta_clusters[meta_clusters$prog_gps =='best',]
      cat('no. of cells in best tumor subtypes=', table(meta_clusters_goodCells$prog_gps), '\n')
      meta_clusters_badCells <- meta_clusters[meta_clusters$prog_gps =='worst',]
      cat('no. of cells in worst tumor subtypes=', table(meta_clusters_badCells$prog_gps), '\n')

      ## ======================
      ## renaming meta-cluster
      ## ======================
      meta_clusters$prog_gps <- ifelse(meta_clusters$prog_gps == 'best', 'GP', 'PP')
      meta_clusters$phenotype <- paste0(meta_clusters$prog_gps,meta_clusters$phenotype)
    }else{
      meta_clusters$phenotype <- paste0('phenotype',meta_clusters$phenotype)
    }

    cat('meta cluster names : ', names(table(meta_clusters$phenotype)), '\n')
    cat('meta cluster size : ', table(meta_clusters$phenotype), '\n')

    ## ======================
    ## cleaning
    ## ======================
    #colnames(meta_clusters)
    meta_clusters <- meta_clusters[, c("cell_ids","core_ids", "phenotype")]

    ## ======================
    ## saving
    ## ======================
    if(prognosis_delineation)save(meta_clusters, file = file.path(output_dir, 'final_meta_clusters_by_prog.RData'))
    if(!prognosis_delineation)save(meta_clusters, file = file.path(output_dir, 'final_meta_clusters.RData'))

    ## ======================
    ## tabulate cluster size in cell counts and fractions
    ## ======================
    freq <- as.matrix(table(meta_clusters$phenotype))
    freq <- data.frame(Phenotype =rownames(freq), cell_counts=freq[,1])
    fractions <- prop.table(table(meta_clusters$phenotype))
    freq$fractions <- fractions
    write.csv(freq, file = file.path(output_dir,'final_meta_cluster_size.csv'), row.names = FALSE)
  }else if(cell_batch_merging == FALSE){
    if(prognosis_delineation)load(file = file.path(output_dir, 'final_meta_clusters_by_prog.RData'))
    if(!prognosis_delineation)load(file = file.path(output_dir, 'final_meta_clusters.RData'))

  }#end cell_batch_merging


}# end function

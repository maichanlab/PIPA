#' Post-PIPA analysis - prognostic association analysis based on cell phenotype composition of tumors
#'
#' Post-PIPA analysis - prognostic association analysis based on cell phenotype composition of tumors
#' Tumor subtyping is performed using \link[Rphenograph]{Rphenograph}, cell phenotype composition of PP,IP&GP as input
#' @param surv_data A data.frame: tumors in rows (tumor_ids in row.names), columns must include: "time", "cens";
#' if there are additional columns, they are treated as confounders
#' @param phenotype_dir A file.path string of input directory storing the PIPA phenotype assignment of all cells of
#' interest
#' @param output_dir A file.path string of output directory
#' @param phenotype_fnm A character string of input file name for cell phenotype data
#' @param density_data A (long format) data.frame containing cell density data (for heatmap viz):
#' samples/tumors in rows with 4 columns named as: "phenotype", "area", "tumor_ids", "density"
#' @param seed An integer of seed for reproducibility
#' @param min_cluster_size An integer indicating the minimum cluster size i.e. cluster/tumor subtypes containing
#' fewer tumors than this number will be excluded in downstream survival analysis
#' i.e. input argument for \code{surv_analysis_fn} \code{Viz_heatmap_bySubtype}
#' @param neighbor_size_byFrac A numerical vector, in fraction, indicating the fraction of total samples/tumors to be
#' used in PhenoGraph
#' @param KM_xlab A character string for labelling the x-axis in Kaplan-Meier curves
#' @param area_split A boolean indicating if the analysis (fraction calculations) will be split by area
#' @import dplyr survival survminer grid RColorBrewer
#' @importFrom tidyr spread
#' @export
composition_surv_analysis <- function(surv_data= NULL,
                                      phenotype_fnm = NULL,
                                      density_data=NULL,
                                      phenotype_dir=NULL, output_dir=NULL,
                                      seed=999, min_cluster_size= 20,
                                      neighbor_size_byFrac= c(0.05, 0.1),
                                      KM_xlab='Survivals',area_split=TRUE){
                                      #, ...){
  phenotype <- tumor_ids <- cluster_no <- median <- NULL
  ## ======================
  ## identify cell phenotype result files
  ## ======================
  setwd(phenotype_dir)
  if(is.null(phenotype_fnm)){
    phenotype_fnm <- 'cell_phenotype.RData'
  }
  ## ======================
  ## identify confounders
  ## ======================
  confounders <- setdiff(colnames(surv_data), c('cens','time'))

  ## ======================
  ## complete cases with survival
  ## ======================
  surv_data <- surv_data[complete.cases(surv_data[,c("cens","time")]),]

  ## ---------------
  ## color palette
  ## ---------------
  col_vec <- c(brewer.pal(9, 'Set1')[-6], brewer.pal(8, 'Set2'), brewer.pal(12,'Set3')[-c(2,12)],
               brewer.pal(8, 'Accent')[-4],brewer.pal(12, 'Paired')[-11])


  ## ======================
  ## load cell phenotype data
  ## ======================
  setwd(phenotype_dir)
  all_cells <- get(load(phenotype_fnm))
  area <- 'overall'
  all_cells$area <- area

  ## ======================
  ## subsetting dens data by complete cases with survival
  ## ======================
  all_cells <- all_cells[all_cells$tumor_ids %in% rownames(surv_data), ]
  cat('#tumors with available survival data=',length(unique(all_cells$tumor)), '\n')
  ## ======================
  ## split by area
  ## ======================
  if(area_split){
    cellids_temp <- sapply(strsplit(x=all_cells$cell_ids, split = 'XY') , '[[', 2)
    possibleError <- tryCatch(
      area <- sapply(strsplit(x=cellids_temp, split = '_') , '[[', 3),
      error=function(e) e
    )
    if(!inherits(possibleError, "error")){
      all_cells$area <- area
    }
  }

  cat('#cells in areas: ',  names(table(all_cells$area)), '\n')
  cat('#cells in areas: ',  table(all_cells$area), '\n')

  for(aa in unique(all_cells$area)){
    cat('At area: ', aa, '\n')
    ## ======================
    ## subsetting for area
    ## ======================
    all_cells_sub <- all_cells[all_cells$area == aa,]

    ## ======================
    ## counting cell phenotypes by cores
    ## ======================
    cell_count_byTumor <- all_cells_sub %>% group_by(phenotype, tumor_ids) %>% count()

    ## ======================
    ## convert cellcounts data to wide format
    ## ======================
    cell_count_byTumor_wide <- spread(cell_count_byTumor[,c('tumor_ids',"phenotype","n")], key='phenotype', value='n')

    ## ======================
    ## combine cell counts with surv data
    ## ======================
    cell_count_byTumor_wide <- merge(cell_count_byTumor_wide, surv_data, by.x = 'tumor_ids', by.y='row.names')


    dens_colids <- grep(x=colnames(cell_count_byTumor_wide),
                        pattern = paste0(unique(all_cells$phenotype), collapse = '|'), value = FALSE)
    for(dd in dens_colids){
      cell_count_byTumor_wide[,dd][is.na(cell_count_byTumor_wide[,dd])] <- 0
    }

    ## ======================
    ## calculate phenotype fractions in individual tumor
    ## ======================
    data <- cell_count_byTumor_wide[, dens_colids, drop=FALSE]
    rownames(data) <- cell_count_byTumor_wide$tumor_ids
    rs <- rowSums(data, na.rm = TRUE)
    fractions <- cbind(data)/ (rs)


    ## ======================
    ## create output subdirectory
    ## ======================
    output_subdir <- file.path(output_dir, paste0(aa,'_compositionSurv'))
    dir.create(output_subdir)

    ## ======================
    ## randomize sample ordering
    ## ======================
    set.seed(seed)
    rand_order <- sample(x = nrow(fractions), size = nrow(fractions))
    fractions <- fractions[rand_order,]

    ## ======================
    ## exclude tumor with equal fractions of all phenotypes: causing NA value during correlation calculation
    # as pearson distance is used in Hierarchical clustering
    ## ======================
    # sd_test <- apply(X = as.data.frame(fractions), MARGIN = 1, sd)
    # equalFrac_ids <- which(sd_test==0)
    # if(length(equalFrac_ids)>0) warning('Tumors with equal fractions of phenotypes (to be excluded):',
    #                                         names(equalFrac_ids))
    # fractions <- fractions[sd_test!=0,]

    ## ======================
    ## Phenograph clustering of tumors: subtyping
    ## ======================
    for(ff in neighbor_size_byFrac){
      neighbor_size <- floor(ff * nrow(fractions))
      set.seed(seed)
      Rphenograph_out <- Rphenograph(fractions, k = neighbor_size)

      tumor_clusters <- data.frame(tumor_ids= rownames(fractions), cluster_no=Rphenograph_out[[2]]$membership)

      output_subdir2 <- file.path(output_subdir,paste0('neigbor',ff))
      dir.create(output_subdir2)
      write.csv(tumor_clusters, file = file.path(output_subdir2,'cluster_k.csv'), row.names = FALSE)

      res_k2 <- merge(tumor_clusters, fractions, by.x='tumor_ids', by.y='row.names')
      write.csv(res_k2 , file=file.path(output_subdir2,'Frac_byTumor.csv'), row.names = FALSE)

      meanFrac_by_cluster <- res_k2 %>% select(-tumor_ids) %>% group_by(cluster_no) %>% summarise_all(mean)
      #medianFrac_by_cluster <- res_k2 %>% select(-tumor_ids) %>% group_by(cluster_no) %>% summarise_all(median)
      medianFrac_by_cluster <- res_k2 %>% select(-tumor_ids) %>% group_by(cluster_no) %>% summarise_all(quantile, probs=0.5)
      write.csv(meanFrac_by_cluster , file=file.path(output_subdir2,'meanFrac_by_cluster.csv'), row.names = FALSE)
      write.csv(medianFrac_by_cluster , file=file.path(output_subdir2,'medianFrac_by_cluster.csv'), row.names = FALSE)


      ## ======================
      ## survival analysis using identified composition-based tumor subtypes
      ## ======================
      surv_analysis_fn(surv_data=surv_data, subtype=tumor_clusters,
                       output_dir=output_subdir2, min_cluster_size=min_cluster_size)

      #### phenotype composition heat-map
      fractions_sub <- fractions[rownames(fractions) %in% tumor_clusters$tumor_ids,]
      fractions_sub$tumor_ids <- rownames(fractions_sub)

      Viz_heatmap_bySubtype(data = fractions_sub, density_data=density_data, subtype=tumor_clusters,
                               output_dir=output_subdir2,feature_order = c('PP','IP','GP'),
                               min_cluster_size=min_cluster_size,
                               value_name='Fraction', row_title='Phenotype')
    }

  }# end all areas
}

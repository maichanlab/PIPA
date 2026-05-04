#' Heat-map visualization: distribution of features grouped by tumor subtypes
#'
#' Heat-map visualization: distribution of features (e.g. cellular features, PIPA-dervied phenotype density)
#' grouped by tumor subtypes with optional annotation of cell density data
#' @param output_dir A file.path to the output directory
#' @param data A data.frame containing tumor/case-level values:
#' tumors/cases in rows with 1+N columns named as: "tumor_ids","value1", "value2" etc for N features
#' @param density_data A (long format) data.frame containing cell density data:
#' tumors/cases in rows with 4 columns named as: "phenotype", "area", "tumor_ids", "density"
#' @param subtype A data.frame containing tumor subtyping data:
#' tumors/cases in rows with 2 columns named as: "tumor_ids","cluster_no"
#' @param min_cluster_size An integer indicating the minimum cluster size (i.e. no. of cases/tumors) for which
#' a tumor subtype containing tumors less than this number will be excluded in the heatmap#'
#' @param log10_dens A boolean: 1 if density is log10 transformed, 0 if otherwise
#' @param value_name A character indicating the name argument passed to \link[ComplexHeatmap]{Heatmap}
#' @param row_title A character indicating the row_title argument passed to \link[ComplexHeatmap]{Heatmap}
#' @param feature_order A vector of characters represent the feature levels (i.e., ordering of plotting)
#' @param plot_width An integer to be passed to the 'width' argument in
#' \link[ggplot2]{ggsave}
#' @param plot_height An integer to be passed to the 'height' argument in \link[ggplot2]{ggsave}
#' @param column_scale A boolean indicating if heatmap values are column/marker-wise scaled
#' @param plot_bnm A character to be appended to the output plot filename
#' @import dplyr grid RColorBrewer
#' @export
Viz_heatmap_bySubtype <- function(data = NULL, density_data=NULL, subtype=NULL,
                                  output_dir=NULL, min_cluster_size= 20, log10_dens=TRUE,
                                  value_name='Value', row_title='Feature',
                                  feature_order=NULL,
                                  plot_bnm=NULL, column_scale=FALSE,
                                  plot_width=12, plot_height=10){

  if(!is.null(plot_bnm)){plot_bnm <- paste0("_",plot_bnm)}
  ## ---------------
  ## color palette
  ## ---------------
  col_vec <- c(brewer.pal(9, 'Set1')[-6], brewer.pal(8, 'Set2'), brewer.pal(12,'Set3')[-c(2,12)],
               brewer.pal(8, 'Accent')[-4],brewer.pal(12, 'Paired')[-11])


  no_cluster <- length(unique(subtype$cluster_no))

  subtype$cluster_no <- factor(subtype$cluster_no)
  no_cluster <- length(unique(subtype$cluster_no))
  col_vec2 <- col_vec[c(1:no_cluster)]
  names(col_vec2) <- as.character(levels(subtype$cluster_no))

  #### exclude cluster with small size
  freq <- data.frame(table(subtype$cluster_no))
  if(max(freq$Freq) < min_cluster_size)stop('ALL clusters are smaller than the minimum required size!\n')
  outlier_clusterids <- freq$Var1[which(freq$Freq < min_cluster_size)]
  subtype <- subtype[!subtype$cluster_no %in% outlier_clusterids,]
  subtype$cluster_no <- droplevels(subtype$cluster_no)
  if( length(unique(subtype$cluster_no))<=1 ){
    warning('Only 1 cluster is found!\n')
    return(NULL)
  }
  rm(freq)


  ## -------------
  ## loop over each area and phenotype
  ## -------------
  if(!is.null(density_data)){
    area <- phenotype <- NULL
    uArea_phenotype <- unique(density_data[,c('area','phenotype')])
  }else{
    uArea_phenotype <- data.frame(area='dummy',phenotype='dummy')
  }

  for(pp in 1:nrow(uArea_phenotype)){
    cat(paste0(uArea_phenotype[pp,], collapse=': '), '\n')

    if(!is.null(density_data)){
      dens_sub <- density_data %>% filter(area == uArea_phenotype[pp,"area"]) %>%
        filter(phenotype == uArea_phenotype[pp, "phenotype"])
      ## find common cases with clustering and density data
      dens_sub <- dens_sub[dens_sub$tumor_ids %in% subtype$tumor_ids,]

      ## exclude cases with no density data
      dens_sub <- dens_sub[complete.cases(dens_sub$density), ]
      data_sub <- data[data$tumor_ids %in% dens_sub$tumor_ids, ]
      tumor_clusters_sub <- subtype[subtype$tumor_ids %in% dens_sub$tumor_ids, ]
    }else{
      data_sub <- data;
      tumor_clusters_sub <- subtype
    }

    ## ordered by clustering assignment
    ## order by clustering number ascendingly
    tumor_clusters_sub <- tumor_clusters_sub[order(tumor_clusters_sub$cluster_no),]
    tumor_clusters_sub$tumor_ids <- as.character(tumor_clusters_sub$tumor_ids)

    ## subsetting for data values
    data_sub <- data_sub[data_sub$tumor_ids %in% tumor_clusters_sub$tumor_ids,]

    ## matching tumor ids ordering of P-value
    data_sub$tumor_ids <- as.character(data_sub$tumor_ids)
    data_sub <- data_sub[match(tumor_clusters_sub$tumor_ids, data_sub$tumor_ids),]
    if(! identical(data_sub$tumor_ids, tumor_clusters_sub$tumor_ids) )
      stop('Tumor ids between data and consensus hierarachical clustering are NOT matched\n')
    ## matching tumor ids ordering of density
    if(!is.null(density_data)){
      dens_sub$tumor_ids <- as.character(dens_sub$tumor_ids)
      dens_sub <- dens_sub[dens_sub$tumor_ids %in% tumor_clusters_sub$tumor_ids,]
      dens_sub <- dens_sub[match(tumor_clusters_sub$tumor_ids, dens_sub$tumor_ids),]
      if(!identical(dens_sub$tumor_ids, tumor_clusters_sub$tumor_ids) )
        warning('Tumor ids between density data and consensus hierarachical clustering are NOT matched\n')

      ## density data check
      if(sum(is.na(dens_sub$density)) > 0){
        invalid_density <- dens_sub$tumor_ids[is.na(dens_sub$density)]
        warning(paste0(invalid_density, collapse = ', '), ': invalid density!!\n')
        break
      }
      if(sum(dens_sub$density > 0) == 0) {
        warning('All tumors have 0 density\n')
        break
      }
      if(log10_dens)dens_sub$density <- log(dens_sub$density +0.1, base = 10)
    }
    tumor_clusters_sub$cluster_no <- as.factor(tumor_clusters_sub$cluster_no)

    ## heat-map column annotation
    ha_column = ComplexHeatmap::HeatmapAnnotation(df = data.frame(Subtype = tumor_clusters_sub$cluster_no),
                                                  col=list(Subtype = col_vec2),
                                                  gap = unit(1.5, "mm"), height = unit(2,"cm"),
                                                  annotation_legend_param = list(Subtype = list(title = "Subtype", title_gp = gpar(fontsize = 18),
                                                                                                labels_gp = gpar(fontsize = 14))))

    if(!is.null(density_data)){
      hbtm = ComplexHeatmap::HeatmapAnnotation(density=ComplexHeatmap::anno_points(dens_sub$density,
                                                                                   height = unit(4, "cm")))
      hbtm@anno_list$density@name <- paste0(uArea_phenotype[pp,"phenotype"],' density in ', uArea_phenotype[pp,"area"])
    }else {hbtm <- NULL}

    ## heat-map
    data_sub[,c('tumor_ids','cluster_no')] <- NULL
    # order features
    if(!is.null(feature_order)){
      feature_order <- feature_order[feature_order%in% colnames(data_sub)]
      data_sub <- data_sub[, match(feature_order,colnames(data_sub))]
    }

    h1 = ComplexHeatmap::Heatmap(t(data_sub), row_title = row_title,
                                 column_title = "Case", name=value_name,
                                 column_dend_reorder = FALSE,
                                 top_annotation = ha_column,cluster_columns=FALSE,
                                 bottom_annotation=hbtm,
                                 show_column_names = FALSE,
                                 heatmap_legend_param = list(title_gp = gpar(fontsize = 16), labels_gp = gpar(fontsize = 14),
                                                             legend_direction = "horizontal",
                                                             legend_height = unit(8, "cm"),legend_width = unit(5, "cm"),
                                                             title_position = "lefttop"))

    if(!is.null(density_data)){
      fileNm <- file.path(output_dir,paste0('heatmap',plot_bnm,"_",paste0(uArea_phenotype[pp,],collapse = '_'),'.pdf'))
      if(log10_dens)fileNm <- gsub(x=fileNm, pattern = '.pdf', replacement = 'log10.pdf')
    }else{
      fileNm <- file.path(output_dir,paste0('heatmap',plot_bnm,'.pdf'))
    }

    pdf(fileNm, width = plot_width, height = plot_height)
    ComplexHeatmap::draw(h1,heatmap_legend_side = "bottom")
    dev.off()


    if(column_scale){
      data_sub_zscore <- scale(data_sub)
      h1 = ComplexHeatmap::Heatmap(t(data_sub_zscore), row_title = row_title,
                                   column_title = "Case", name=value_name,
                                   column_dend_reorder = FALSE,
                                   top_annotation = ha_column,cluster_columns=FALSE,
                                   bottom_annotation=hbtm,
                                   show_column_names = FALSE,
                                   heatmap_legend_param = list(title_gp = gpar(fontsize = 16), labels_gp = gpar(fontsize = 14),
                                                               legend_direction = "horizontal",
                                                               legend_height = unit(8, "cm"),legend_width = unit(5, "cm"),
                                                               title_position = "lefttop"))

      if(!is.null(density_data)){
        fileNm <- file.path(output_dir,paste0('heatmap',plot_bnm,paste0(uArea_phenotype[pp,],collapse = '_'),'_ct_zscore.pdf'))
        if(log10_dens)fileNm <- gsub(x=fileNm, pattern = '.pdf', replacement = 'log10.pdf')
      }else{
        fileNm <- file.path(output_dir,paste0('heatmap',plot_bnm,'_ct_zscore.pdf'))
      }

      pdf(fileNm, width = plot_width, height = plot_height)
      ComplexHeatmap::draw(h1,heatmap_legend_side = "bottom")
      dev.off()
    }## end heat-map with celltype-wise zscores
  }

}

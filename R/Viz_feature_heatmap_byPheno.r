#' Heat-map visualization: distribution of cell features grouped by PIPA-derived phenotypes
#'
#' Heat-map visualization: distribution of cell features grouped by PIPA-derived phenotypes
#' using (1) case/tumor-median feature values, and (2) cellular feature values of
#' PIPA-selected features
#' @param cell_data A data.frame contains cells in rows and features in columns; MUST contain a column named 'cell_ids'
#' as cell identifiers (core_ids+coordinates) and a column named 'tumor_ids' for tumor/case identifiers
#' @param pheno_levels A vector of characters represent the phenotype levels (i.e., ordering of plotting)
#' @param pheno_dir A file.path string of input directory storing the PIPA phenotype assignment of all cells of
#' interest
#' @param plot_width An integer to be passed to the 'width' argument in
#' \link[ggplot2]{ggsave}
#' @param plot_height An integer to be passed to the 'height' argument in
#' \link[ggplot2]{ggsave}
#' @param output_dir A file.path to the output directory
#' @param phenotype_fnm A character string of input file name for cell phenotype data
#' @param column_scale A boolean indicating if heatmap values are column/marker-wise scaled
#' @param filetype A character string indicating the output plot file format, pdf or png
#' @return PDF & PNG files containing the heat maps
#' @import ComplexHeatmap
#' @importFrom plyr mapvalues
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices png
#' @export
Viz_feature_heatmap_byPheno <- function(cell_data= NULL,
                                        phenotype_fnm = "cell_phenotype.RData",
                                        pheno_dir=NULL, output_dir=NULL,
                                        pheno_levels=c('PP','IP','GP'),
                                        column_scale=TRUE,filetype='pdf',
                                        plot_width=NULL, plot_height=8){

  tumor_ids <- phenotype <- median <- NULL

  ## ======================
  ## identify plotting features
  ## ======================
  plot_features <- colnames(cell_data)
  plot_features <- setdiff(plot_features, c('tumor_ids','cell_ids'))
  ## ======================
  ## create output sub-folder
  ## ======================
  output_subdir <- file.path(output_dir,
                             paste0('Viz_heatmap'))
  dir.create(output_subdir)


  ## ======================
  ## load cell phenotype data
  ## ======================
  setwd(pheno_dir)
  all_cells <- get(load(phenotype_fnm))
  if(sum(duplicated(all_cells$cell_ids))>0){
    all_cells <- all_cells[-which(duplicated(all_cells$cell_ids)),]
  }
  ## ======================
  ## append phenotype assignment to cell feature data
  ## ======================
  cell_data <- merge(cell_data, all_cells[,c("cell_ids","phenotype")], by='cell_ids')

  ## ======================
  ## subsetting for features
  ## ======================
  cell_data_features <- cell_data[, grep(x=colnames(cell_data),
                                         pattern = paste0(c(plot_features,'cell_ids','tumor_ids',"phenotype"), collapse = '|'))]

  ## ======================
  ## set level for phenotypes
  ## ======================
  if(is.null(pheno_levels)){
    phenotype_map <- data.frame(phenotype=unique(cell_data_features$phenotype))
    phenotype_map$type <- ifelse(grepl(phenotype_map$phenotype, pattern = 'GP'), 'GP',
                                 ifelse(grepl(phenotype_map$phenotype, pattern = 'PP'), 'PP', 'others'))
    phenotype_map$number <- gsub(x=phenotype_map$phenotype, pattern = 'GP|PP', replacement = '')
    phenotype_map$number <- as.integer(phenotype_map$number)
    phenotype_map$type <- factor(phenotype_map$type, levels=c('PP','GP'))
    phenotype_map <- phenotype_map[order(phenotype_map$type, phenotype_map$number),]
    cell_data_features$phenotype <- factor(cell_data_features$phenotype,
                                           levels = phenotype_map$phenotype)

  }else{
    cell_data_features$phenotype <- factor(cell_data_features$phenotype,
                                           levels = pheno_levels)
  }
  ## ======================
  ## cell level expression
  ## ======================
  cell_data_features_cellLvl <- cell_data_features
  rownames(cell_data_features_cellLvl) <- cell_data_features_cellLvl$cell_ids
  cell_data_features_cellLvl$tumor_ids <- NULL
  cell_data_features_cellLvl$cell_ids <- NULL

  ## ======================
  ## compute tumor median expression
  ## ======================
  cell_data_features_tumorLvl <- cell_data_features
  cell_data_features_tumorLvl$cell_ids <- NULL
  # cell_data_features_tumorLvl <- cell_data_features_tumorLvl %>% group_by(tumor_ids,phenotype) %>%
  #   summarise_all(median, na.rm=TRUE)
  cell_data_features_tumorLvl <- cell_data_features_tumorLvl %>% group_by(tumor_ids,phenotype) %>%
    summarise_all(quantile,probs=0.5, na.rm=TRUE)
  cell_data_features_tumorLvl$tumor_ids <- NULL

  ## ======================
  ## sort by phenotypes
  ## ======================
  cell_data_features_tumorLvl <- cell_data_features_tumorLvl[order(cell_data_features_tumorLvl$phenotype),]
  cell_data_features_cellLvl <- cell_data_features_cellLvl[order(cell_data_features_cellLvl$phenotype),]

  ## ======================
  ## color annotation for phenotypes
  ## ======================
  if(length(unique(cell_data_features$phenotype))==3){
    col_vec <- c('#F8766D','#00BA38','#619CFF')
  }else{
    col_vec <- c(brewer.pal(9, 'Set1'), brewer.pal(8, 'Set2'), brewer.pal(12,'Set3'),
                 brewer.pal(8, 'Accent'),brewer.pal(12, 'Paired'))
  }
  col_vec2 <- col_vec[1:length(unique(cell_data_features$phenotype))]
  names(col_vec2) <- unique(cell_data_features$phenotype)[order(unique(cell_data_features$phenotype))]

  ########################################
  ## plotting: all cells
  ########################################
  pheno_anno = HeatmapAnnotation(df = data.frame(Phenotype = cell_data_features_cellLvl$phenotype),
                                 col=list(Phenotype = col_vec2),
                                 annotation_legend_param = list(Phenotype = list(title = "Phenotype",
                                                                                 title_gp = gpar(fontsize = 16),
                                                                                 labels_gp = gpar(fontsize = 16),
                                                                                 nrow = 1, by_row = TRUE)))

  ## order features
  temp_cellLvl <- cell_data_features_cellLvl[,match(plot_features,colnames(cell_data_features_cellLvl)), drop=FALSE]

  name="Cellular value (raw)"

  if(column_scale){
    temp_cellLvl <- scale(temp_cellLvl)
    name="Cellular value (scaled)"
  }

  ## heat map
  ht1 = Heatmap(matrix = as.matrix(t(temp_cellLvl)),
                name = name, cluster_rows = FALSE,
                cluster_columns = FALSE,
                row_title = 'Features',column_title= 'Cells',show_column_names=FALSE,
                row_title_gp = gpar(fontsize = 20),row_names_gp = gpar(fontsize = 20),
                column_title_gp = gpar(fontsize = 20),
                top_annotation = pheno_anno,
                heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(5, "cm"),
                                            title_gp = gpar(fontsize = 16), title_position ='topleft',
                                            labels_gp = gpar(fontsize = 16)))

  fnm <- file.path(output_subdir, 'cellLvl.pdf')
  if(grepl(filetype, pattern = 'png', ignore.case = T))fnm <- file.path(output_subdir, 'cellLvl.png')


  if(is.null(plot_width)){
    plot_width <- 8
    if(nrow(temp_cellLvl) > 300000)plot_width <- 12
    if(nrow(temp_cellLvl) > 600000)plot_width <- 14
  }
  if(grepl(filetype, pattern = 'png', ignore.case = T)){
    png(fnm)
    draw(ht1, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
    dev.off()
  }else{
    pdf(file = file.path(fnm))
    draw(ht1, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
    dev.off()
  }

  ########################################
  ## plotting: all tumors
  ########################################
  pheno_anno = HeatmapAnnotation(df = data.frame(Phenotype = cell_data_features_tumorLvl$phenotype),
                                 col=list(Phenotype = col_vec2),
                                 annotation_legend_param = list(Phenotype = list(title = "Phenotype",
                                                                                 title_gp = gpar(fontsize = 16),
                                                                                 labels_gp = gpar(fontsize = 16),
                                                                                 nrow = 1, by_row = TRUE)))
  ## order features
  temp_tumorLvl <- cell_data_features_tumorLvl[,match(plot_features,colnames(cell_data_features_tumorLvl)), drop=FALSE]
  name="Case median value (raw)"

  if(column_scale){
    temp_tumorLvl <- scale(temp_tumorLvl)
    name="Case median value (scaled)"
  }

  ht1 = Heatmap(matrix = as.matrix(t(temp_tumorLvl)),
                name = name, cluster_rows = FALSE,
                cluster_columns = FALSE,
                row_title = 'Features', column_title = 'Cases',show_column_names=FALSE,
                row_title_gp = gpar(fontsize = 20),row_names_gp = gpar(fontsize = 20),
                column_title_gp = gpar(fontsize = 20),
                top_annotation = pheno_anno,
                heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(5, "cm"),
                                            title_gp = gpar(fontsize = 16),  title_position ='topleft',
                                            labels_gp = gpar(fontsize = 16)))
  if(is.null(plot_width)){
    plot_width <- 8
    if(nrow(temp_tumorLvl) > 300000)plot_width <- 12
    if(nrow(temp_tumorLvl) > 600000)plot_width <- 14
  }

  fnm <- file.path(output_subdir, 'tumorLvl.pdf')
  pdf(fnm, width = plot_width, height = plot_height)
  draw(ht1, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
  dev.off()



  ########################################
  ## plotting: PP vs GP cells
  ########################################
  cell_data_features_PP <- cell_data_features_GP <- NULL
  if( sum(grep(x=cell_data_features_cellLvl$phenotype,  pattern = 'GP|PP')) >0 ){
    cell_data_features_GP <- cell_data_features_cellLvl %>% filter(grepl(x = phenotype, pattern = 'GP'))
    print(table(cell_data_features_GP$phenotype))
    cell_data_features_PP <- cell_data_features_cellLvl %>% filter(grepl(x = phenotype, pattern = 'PP'))
    print(table(cell_data_features_PP$phenotype))

    ## ======================
    ## GP cells
    ## ======================
    if(nrow(cell_data_features_GP) >0){
      pheno_anno = HeatmapAnnotation(df = data.frame(Phenotype = cell_data_features_GP$phenotype),
                                     col=list(Phenotype = col_vec2),
                                     annotation_legend_param = list(Phenotype = list(title = "Phenotype",
                                                                                     title_gp = gpar(fontsize = 16),
                                                                                     labels_gp = gpar(fontsize = 16),
                                                                                     nrow = 1, by_row = TRUE)))

      temp_GP <- cell_data_features_GP[, grep(x=colnames(cell_data_features_GP),
                                                  pattern = paste0(plot_features, collapse = '|')), drop=FALSE]
      name="Cellular value (raw)"

      if(column_scale){
        temp_GP <- scale(temp_GP)
        name="Cellular value (scaled)"
      }

      ## order features
      temp_GP <- temp_GP[,match(plot_features,colnames(temp_GP)), drop=FALSE]
      ht2 = Heatmap(matrix = as.matrix(t(temp_GP)),
                    name = name, cluster_rows = FALSE,
                    cluster_columns = FALSE,
                    row_title = 'Features', column_title = 'GP cells',show_column_names=FALSE,
                    row_title_gp = gpar(fontsize = 20),row_names_gp = gpar(fontsize = 20),
                    top_annotation = pheno_anno,
                    heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(5, "cm"),
                                                title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 14)))

      if(is.null(plot_width)){
        plot_width <- 8
        if(nrow(temp_GP) > 300000)plot_width <- 12
        if(nrow(temp_GP) > 600000)plot_width <- 14
      }


      fnm <- file.path(output_subdir, 'cellLvl_GP.pdf')
      if(grepl(filetype, pattern = 'png', ignore.case = T))fnm <- file.path(output_subdir, 'cellLvl_GP.png')


      if(grepl(filetype, pattern = 'png', ignore.case = T)){
        png(fnm)
        draw(ht2, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
        dev.off()
      }else{
        pdf(file = file.path(fnm))
        draw(ht2, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
        dev.off()
      }
    }


    ## ======================
    ## PP cells
    ## ======================
    if(nrow(cell_data_features_PP) >0){
      pheno_anno = HeatmapAnnotation(df = data.frame(Phenotype = cell_data_features_PP$phenotype),
                                     col=list(Phenotype = col_vec2),
                                     annotation_legend_param = list(Phenotype = list(title = "Phenotype",
                                                                                     title_gp = gpar(fontsize = 16),
                                                                                     labels_gp = gpar(fontsize = 16),
                                                                                     nrow = 1, by_row = TRUE)))

      temp_PP <- cell_data_features_PP[, grep(x=colnames(cell_data_features_PP),
                                                pattern = paste0(plot_features, collapse = '|')), drop=FALSE]
      if(column_scale)temp_PP <- scale(temp_PP)

      ## order features
      temp_PP <- temp_PP[,match(plot_features,colnames(temp_PP)), drop=FALSE]
      colnames(temp_PP)
      name="Cellular value (raw)"

      if(column_scale){
        temp_PP <- scale(temp_PP)
        name="Cellular value (scaled)"
      }

      ht2 = Heatmap(matrix = as.matrix(t(temp_PP)),
                    name = name, cluster_rows = FALSE,
                    cluster_columns = FALSE,
                    row_title = 'Features', column_title = 'PP cells',show_column_names=FALSE,
                    row_title_gp = gpar(fontsize = 20),row_names_gp = gpar(fontsize = 20),
                    top_annotation = pheno_anno,
                    heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(5, "cm"),
                                                title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 14)))

      plot_width <- 8
      if(nrow(temp_PP) > 300000)plot_width <- 12
      if(nrow(temp_PP) > 600000)plot_width <- 14

      fnm <- file.path(output_subdir, 'cellLvl_PP.pdf')
      if(grepl(filetype, pattern = 'png', ignore.case = T))fnm <- file.path(output_subdir, 'cellLvl_PP.png')


      if(grepl(filetype, pattern = 'png', ignore.case = T)){
        png(fnm)
        draw(ht2, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
        dev.off()
      }else{
        pdf(file = file.path(fnm))
        draw(ht2, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
        dev.off()
      }
    }


  }#end PP & GP heatmaps


}# end of function

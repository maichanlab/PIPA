#' Box-plot visualization: distribution of cell features grouped by PIPA-derived phenotypes
#'
#' Box-plot visualization: distribution of cell features grouped by PIPA-derived phenotypes
#' @param cell_data A data.frame contains cells in rows and features in columns; MUST contain a
#' column named 'cell_ids' as cell identifiers (core_ids & XY-coordinates) and other columns
#' are treated as features for plotting
#' @param pheno_dir A file.path string of input directory storing the PIPA phenotype assignment
#' of all cells of interest
#' @param output_dir A file.path string of output directory
#' @param phenotype_fnm A character string of input file name for cell phenotype data
#' @param pheno_levels A vector of characters represent the phenotype levels (i.e., ordering of plotting)
#' @param outlier_removal A boolean indicating if outliers (outside IQR) are included
#' @param plot_bnm A character string identifier appended to the output folder name
#' @param ncol An integer indicating the number of columns (ncol argument) in
#' \link[ggplot2]{facet_wrap}
#' @param plot_width An integer to be passed to the 'width' argument in
#' \link[ggplot2]{facet_wrap}
#' @param plot_height An integer to be passed to the 'height' argument in
#' \link[ggplot2]{ggsave}
#' @return A PDF file containing the box-plot
#' @importFrom RColorBrewer brewer.pal
#' @import dplyr
#' @importFrom tidyr gather
#' @export
Viz_feature_boxplot_byPheno <- function(cell_data= NULL,
                                            phenotype_fnm = NULL,
                                            pheno_dir=NULL, output_dir=NULL,
                                            plot_bnm='PIPA_features',
                                            pheno_levels =c('PP','IP','GP'),
                                            outlier_removal=TRUE, ncol=1,
                                            plot_width=8, plot_height=8){
  key <- value <- all_of <- phenotype <- tumor_ids <- NULL
  ## ======================
  ## identify cell phenotype result files
  ## ======================
  setwd(pheno_dir)
  if(is.null(phenotype_fnm)){
    phenotype_fnm <- 'cell_phenotype.RData'
  }

  ## ======================
  ## identify plotting features
  ## ======================
  plot_features <- colnames(cell_data)
  plot_features <- setdiff(plot_features, c('cell_ids','tumor_ids'))
  ## ======================
  ## create output sub-folder
  ## ======================
  output_subdir <- file.path(output_dir,paste0('Viz_boxplot'))
  dir.create(output_subdir)


  ## ======================
  ## load cell phenotype data
  ## ======================
  setwd(pheno_dir)
  phenotype_data <- get(load(phenotype_fnm))

  ## ======================
  ## set pheno_levels
  ## ======================
  if(is.null(pheno_levels)){
    pheno_levels <- data.frame(phenotype=unique(phenotype_data$phenotype), stringsAsFactors = FALSE)
    pheno_levels$pheno_ids <- sapply(strsplit(x=pheno_levels$phenotype, 'PP|GP'),'[[', 2)
    pheno_levels$pheno_ids<- as.integer(pheno_levels$pheno_ids )
    pheno_levels$pheno_type <- ifelse(grepl(x=pheno_levels$phenotype, pattern = 'PP'), 'PP',
                                  ifelse(grepl(x=pheno_levels$phenotype, pattern = 'GP'), 'GP','others'))
    pheno_levels <- pheno_levels[order(pheno_levels$pheno_type, pheno_levels$pheno_ids), ]
    pheno_levels <- pheno_levels$phenotype
  }
  ## ======================
  ## append phenotype assignment to cell feature data
  ## ======================
  cell_data <- merge(cell_data, phenotype_data[,c("cell_ids","phenotype")], by='cell_ids')

  ## ======================
  ## subsetting for features
  ## ======================
  cell_data_features <- cell_data[, grep(x=colnames(cell_data),
                                        pattern = paste0(c(plot_features,"phenotype"), collapse = '|'))]



  ## ======================
  ## convert to long format
  ## ======================
  cell_feature_long <- cell_data_features %>% gather(key,value,all_of(plot_features))

  ## ======================
  ## exlude outliers
  ## ======================
  if(outlier_removal){
    filtered_data <- c()
    for(ff in plot_features){
      temp<- cell_feature_long[cell_feature_long$key == ff, ]
      IQR <-IQR(temp$value, na.rm = TRUE)
      Q1 <- quantile(temp$value, prob=0.25, na.rm = TRUE)
      Q3 <- quantile(temp$value, prob=0.75, na.rm = TRUE)
      cell_feature_long2 <- cell_feature_long %>% filter(key==ff) %>%
        filter(value < (Q3+1.5*IQR)) %>% filter(value > (Q1-1.5*IQR))
      filtered_data <- rbind(filtered_data, cell_feature_long2)
    }# end plot_features
    cell_feature_long <- filtered_data
  }# end outlier_removal
  ## ======================
  ## set level for phenotypes
  ## ======================
  cell_feature_long$phenotype <- factor(cell_feature_long$phenotype,
                                             levels = pheno_levels)
  cell_feature_long$phenotype <- droplevels(cell_feature_long$phenotype)
  ## ======================
  ## color annotation for phenotypes
  ## ======================
  if(length(unique(cell_feature_long$phenotype))==3){
    col_vec <- c('#F8766D','#00BA38','#619CFF')
  }else{
    col_vec <- c(brewer.pal(9, 'Set1'), brewer.pal(8, 'Set2'), brewer.pal(12,'Set3'),
                 brewer.pal(8, 'Accent'),brewer.pal(12, 'Paired'))
  }
  col_vec2 <-  col_vec[1:length(unique(cell_feature_long$phenotype))]
  names(col_vec2) <- levels(cell_feature_long$phenotype)

  ## ======================
  ## box-plot
  ## ======================
  pl=ggplot(data=cell_feature_long, aes(x=phenotype, y=value, fill=phenotype))+
    scale_fill_manual(values = col_vec2)+
    geom_boxplot(outlier.color = 'NA')+
    facet_wrap(key ~., scales = 'free_y', ncol = ncol)+
    theme_bw()+
    theme(axis.text = element_text(size = 26, face='bold'),
          axis.title = element_text(size=26, face='bold'),
          axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1),
          strip.text = element_text(size = 28, face='bold'),
          legend.position = 'none')+
    ylab('Values')+xlab('Phenotype')


  fnm <- file.path(output_subdir, paste0(plot_bnm,'.pdf'))
  if(outlier_removal) fnm <- gsub(x=fnm, pattern = '.pdf', replacement = '_rm.outliers.pdf')
  ggsave(plot = pl, filename = fnm, width = plot_width, height = plot_height)


}# end of function

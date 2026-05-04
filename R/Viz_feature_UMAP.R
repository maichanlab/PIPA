#' UMAP visualization: relatedness of cells in the PIPA-selected feature space
#'
#' UMAP visualization: relatedness of cells in the PIPA-selected feature space - annotation by
#' (1) PIPA-derived phenotypes, and (2) cellular feature values of PIPA-selected features
#' @param cell_data A data.frame contains cells in rows and features in columns; MUST contain a column named 'cell_ids'
#' as cell identifiers (core_ids+coordinates) and all other columns will be used as cell features for UMAP dimensional
#' reduction
#' @param downspl_size An integer indicating the number of randomly selected cells
#' (downspl_size cells from PP, IP and GP individually) to be used for UMAP dimensional reduction
#' @param pheno_dir A file.path string of input directory storing the PIPA phenotype assignment of all cells of
#' interest
#' @param seed An integer used as the seed for generating the downsampled cells
#' @param log10_featureNms A vector of characters of feature names whose values to be log10-transformed
#' prior to scaling (z-score); same setting as \code{meta_cluster} should be used
#' @param plot_width An integer to be passed to the 'width' argument in
#' \link[ggplot2]{ggsave}
#' @param plot_height An integer to be passed to the 'height' argument in
#' \link[ggplot2]{ggsave}
#' @param output_dir A file.path to the output directory
#' @param phenotype_fnm A character string of input file name for cell phenotype data
#' @return PDF files containing the UMAP plots
#' @importFrom plyr mapvalues
#' @import umap
#' @export
Viz_feature_UMAP <- function(cell_data= NULL,downspl_size=1000,
                             phenotype_fnm = "cell_phenotype.RData",
                             seed=123,log10_featureNms=NULL,
                             pheno_dir=NULL, output_dir=NULL,
                             plot_width=8, plot_height=8){
  umap_dim1 = umap_dim2=phenotype=NULL
  ## ======================
  ## specify color scheme
  ## ======================
  col_vec <- c('#F8766D','#00BA38','#619CFF')
  col_map <- data.frame(phenotype=c('PP','IP','GP'), stringsAsFactors = TRUE)
  col_map$rowids <- 1:nrow(col_map)
  col_map$color <- col_vec[1:nrow(col_map)]
  col_map$phenotype <- factor(col_map$phenotype, levels = c('PP','IP','GP'))
  col_map <- col_map[order(col_map$phenotype),]

  ## ======================
  ## create output sub-folder
  ## ======================
  output_subdir <- file.path(output_dir,
                             paste0('Viz_UMAP'))
  dir.create(output_subdir)

  ## ======================
  ## load cell phenotype data
  ## ======================
  setwd(pheno_dir)
  pheno_data <- get(load(phenotype_fnm))
  if(sum(duplicated(pheno_data$cell_ids))>0){
    pheno_data <- pheno_data[-which(duplicated(pheno_data$cell_ids)),]
  }
  ## ======================
  ## append phenotype assignment to cell feature data
  ## ======================
  m <- merge(cell_data, pheno_data[,c("cell_ids","phenotype")], by='cell_ids')

  ## ======================
  ## count no. of cells in each phenotype
  ## ======================
  freq <- as.data.frame(as.matrix(table(m$phenotype)))
  colnames(freq) <- 'cell_counts'
  freq$pheno <- rownames(freq)

  ########################################
  ## downsampling
  ########################################
  m_IP <- m[m$phenotype=='IP',]
  if(nrow(m_IP) > 0){
    set.seed(seed)
    down_IP <- sample(nrow(m_IP), downspl_size, replace = FALSE)
    m_IP <- m_IP[down_IP, ]
  }
  m_GP <- m[m$phenotype=='GP',]
  if(nrow(m_GP) > 0){
    set.seed(seed)
    down_GP <- sample(nrow(m_GP), downspl_size, replace = FALSE)
    m_GP <- m_GP[down_GP, ]
  }
  m_PP <- m[m$phenotype=='PP',]
  if(nrow(m_PP) > 0){
    set.seed(seed)
    down_PP <- sample(nrow(m_PP), min(freq$cell_counts[freq$pheno=='PP'],downspl_size), replace = FALSE)
    m_PP <- m_PP[down_PP, ]
  }

  m <- rbind(m_IP, m_GP, m_PP)
  print(table(m$phenotype))

  ## ======================
  ## extract cell feature values
  ## ======================
  selected_features <- setdiff(colnames(cell_data), c('cell_ids','tumor_ids'))
  m2 <- as.data.frame(m)[,colnames(m) %in% selected_features]
  rownames(m2) <- m$cell_ids

  ## ======================
  ## data transformation
  ## ======================
  for (cc in 1:ncol(m2)){
    feature_cc <- colnames(m2)[cc]
    if(feature_cc %in% log10_featureNms){
      cat('LOG10:', feature_cc, '\n')
      m2[,cc] <-log10_transform(m2[,cc], method = 'log10')
    }
    m2[,cc] <- scale(m2[,cc])
  }
  ## ======================
  ## UMAP analysis
  ## ======================
  set.seed(seed)
  umap_res <- umap(d=m2)

  ## ======================
  ## plotting with phenotype annotation
  ## ======================
  ## check phenotype data ordered
  umap_res_df <- as.data.frame(umap_res$layout)
  stopifnot(identical(rownames(umap_res_df), m$cell_ids))

  ## append phenotyp data to umap results
  umap_res_df$phenotype <- m$phenotype
  umap_res_df$phenotype <- factor(umap_res_df$phenotype, levels = col_map$phenotype)
  colnames(umap_res_df) <- mapvalues(x =colnames(umap_res_df), from = c('V1','V2'),
                                     to = c('umap_dim1', 'umap_dim2'))

  save(umap_res_df, file=file.path(output_subdir, 'umap_results.RData'))

  ## plotting
  umap_plot <-ggplot(data = umap_res_df, aes(x=umap_dim1, y = umap_dim2, color=phenotype))+
    scale_color_manual(values = c(col_map$color))+
    geom_point(size=0.4)+theme_bw()+
    theme(legend.position = 'bottom',
          #axis.title = element_text(size=7, family = 'Arial'),
          #axis.text = element_text(size=7, family = 'Arial'),
          axis.title = element_text(size=7),
          axis.text = element_text(size=7),
          axis.ticks =element_blank(),
          panel.border = element_rect(size = 0.3),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    xlab('uMAP dim1') + ylab('uMAP dim2')
  ggsave(plot=umap_plot, filename = file.path(output_subdir, 'umap_byPheno.pdf'))

  ## ======================
  ## plotting with cell feature annotation
  ## ======================
  ## assign cells into decile categories for better viz
  m3 <- m2
  for (cc in 1:ncol(m3)){
    feature_cc <- colnames(m3)[cc]
    temp <- m3[,cc]
    q <- quantile(temp, probs = seq(0, 1, by = 0.1), na.rm = TRUE)
    q <- unique(q)
    cat <- cut(temp,  breaks=q, include.lowest=TRUE)
    cat(feature_cc, '\n')
    cat(table(cat),'\n')
    m3[,cc] <- as.integer(cat)
    cat(table(m3[,cc]),'\n')
  }
  m3 <- m3
  m3$cell_ids <- rownames(m3)

  ## merge decile-categorized feature value with UMAP results
  umap_res_df$cell_ids <- rownames(umap_res_df)
  umap_res_df <- left_join(umap_res_df, m3, by='cell_ids')

  ## plotting
  for(ff in setdiff(colnames(m3),'cell_ids')){
    col_ff <- scale_colour_gradientn(colours = c('red','green'), limits=c(1, 10))
    target_plot <- ggplot(data = umap_res_df, aes_string(x='umap_dim1', y = 'umap_dim2', color=ff))+
      geom_point(size=0.8)+ theme_bw()+
      theme(legend.position = 'bottom',
            #title = element_text(size=5.8, family = 'Arial'),
            title = element_text(size=5.8),
            #axis.title = element_text(size=7, family = 'Arial'),
            axis.title = element_text(size=7),
            axis.ticks =element_blank(),
            #axis.text = element_text(size=7, family = 'Arial'),
            axis.text = element_text(size=7),
            panel.border = element_rect(size = 0.3),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())+col_ff+
      xlab('uMAP dim1') + ylab('uMAP dim2')+
      ggtitle(ff)
    ggsave(plot=target_plot, filename = file.path(output_subdir, paste0('umap_by',ff,'_decileColors.pdf')))
  }# end all feature annotation

}# end of function

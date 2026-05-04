#' Optimize parameters
#'
#' Compare and optimize parameters, based on the distribution of (a) %cell, and
#' (b) prevalence of patient/tumor for  PP, IP, and GP phenotypes
#' @param root_dir A file.path pointing to the parent directory of the folders containing
#' results using different param settings
#' @param pheno_dir_path A character string of the directory path to the (final) phenotype & density data
#' @param param_bnm A character string of the key name for identifying the folder names of
#' the experimental setting of interest
#' @param output_dir A file.path string of output directory
#' @param pheno_dens_fnm A character string of input file name for cell phenotype density data
#' (after merging i.e. PP, IP, GP)
#' @param pheno_count_fnm A character string of input file name for cell phenotype count data
#' (after merging i.e. PP, IP, GP)
#' @param tissue_area_list A character vectors indicating the tissue areas of interest for plotting
#' @param plot_bnm A character string identifier to be appended to the output PDF
#' @param plot_width An integer to be passed to the 'width' argument in
#' \link[ggplot2]{facet_wrap}
#' @param plot_height An integer to be passed to the 'height' argument in
#' \link[ggplot2]{ggsave}
#' @return Two PDF files
#' @import dplyr ggplot2
#' @importFrom stats reshape
#' @export
Optimize_Params <-
  function(root_dir=NULL,
           pheno_dir_path=NULL,
           output_dir= NULL,param_bnm = NULL,
           pheno_dens_fnm = 'cell_phenotype_density.RData',
           pheno_count_fnm = 'cell_phenotype.RData',
           tissue_area_list=c('T','S'),
           plot_bnm='',plot_width=8, plot_height=8){

    key <- phenotype <- density <- nonzero <- density_check <- detection_rate <- no_of_cells <- NULL

    ## ======================
    ## color annotation for phenotypes
    ## ======================
    col_vec <- c('PP'='#F8766D','IP'='#00BA38','GP'='#619CFF')


    folder_nms <- list.files(path = root_dir,pattern = param_bnm)
    for(tissue_area in tissue_area_list){
      ########################################
      #### bar plots: percentage of cases containing
      # nonzero cells of PP, IP, and GP phenotypes
      ########################################

      ## ================
      ## load and combined cell density data from all input folders
      ## ================
      dens_all_cutoff <- c()
      for (ff in 1:length(folder_nms)){
        if(!file.exists((file.path(root_dir,folder_nms[ff],pheno_dir_path,pheno_dens_fnm))))next
        setwd(file.path(root_dir,folder_nms[ff],pheno_dir_path))
        temp <- get(load(file=pheno_dens_fnm))
        temp <- temp[temp$area %in% tissue_area,]
        temp$key <- folder_nms[ff]
        #cat(nrow(temp), '\n')
        dens_all_cutoff <- rbind(dens_all_cutoff, temp)
      }

      ## ================
      ## count no. of cases/tumors containing zero & nonzero phenotypes
      ## ================
      case_count_by_phenotypes <- dens_all_cutoff  %>% ungroup() %>%
        group_by(key, phenotype) %>% count(density>0, name = 'nonzero')
      colnames(case_count_by_phenotypes)[grep(x=colnames(case_count_by_phenotypes), pattern = '^density')] <- 'density_check'


      ## ================
      ## compute % cases containing nonzero phenotype density
      ## ================
      perc_nonzero <- case_count_by_phenotypes %>% group_by(key, phenotype) %>%
        summarise(perc_nonzero= nonzero[density_check]/sum(nonzero))
      ## ================
      ## fill in empty phenotype class/cat
      ## ================
      #missing_cat <- data.frame(key='detection_rate0', phenotype='PP', perc_nonzero=0, stringsAsFactors = FALSE)
      #perc_nonzero2 <- rbind(as.data.frame(perc_nonzero), missing_cat)
      ## ================
      ## stacked bar-plot
      ## ================
      perc_nonzero2<-perc_nonzero
      perc_nonzero2$detection_rate <- gsub(x=perc_nonzero2$key, pattern = param_bnm, replacement = '')
      perc_nonzero2$detection_rate <- as.numeric(as.character(perc_nonzero2$detection_rate))
      perc_nonzero2$detection_rate <- as.factor(perc_nonzero2$detection_rate)
      perc_nonzero2$phenotype <- factor(perc_nonzero2$phenotype,
                                        levels = c('PP','IP','GP'))

      ylab <- expression(paste("%cases > 0 cell"))
      pl1 <- ggplot(perc_nonzero2, aes(fill=phenotype, y=perc_nonzero, x=detection_rate)) +
        geom_bar(stat="identity",position=position_dodge())+
        theme_bw()+
        theme(axis.text = element_text(color='black',size = 11),
              axis.title = element_text(color='black',size=11),
              panel.border = element_rect(color='black',size=0.5),
              panel.grid = element_blank(),
              legend.spacing.x = unit(0.2, 'cm'),
              legend.title =  element_text(size = 11),
              legend.text = element_text(size = 11),
              legend.position = 'bottom',
              legend.direction = 'horizontal')+
        ylab(ylab)+ xlab(param_bnm)+
        guides(fill = guide_legend(title = 'phenotype', title.position = "left"))+
        guides(fill = guide_legend(override.aes = list(size = 3)))

      setwd(output_dir)
      ggsave(plot = pl1, width = plot_width, height = plot_height,
             filename = paste0(paste0('prevalence_by',param_bnm, '_in_', tissue_area,'.pdf')))

      ## ==============
      # saving prevalent data
      ## ==============
      save(perc_nonzero2,
           file = file.path(output_dir,paste0('prevalence_by',param_bnm, '_in_', tissue_area,'.RData')))

      ## ==============
      ## Identify the optimal setting
      ## ==============
      perc_nonzero2$detection_rate  <- NULL
      perc_nonzero2 <- as.data.frame(perc_nonzero2)
      dt_wide_prev <- reshape(perc_nonzero2, idvar = c("key"),
                              timevar = "phenotype", direction = "wide")

      ## replace NAs with Zeros
      dt_wide_prev[, names(dt_wide_prev) != "key"] <- lapply(dt_wide_prev[, names(dt_wide_prev) != "key"], function(x) ifelse(is.na(x), 0, x))
      ## clean column names
      colnames(dt_wide_prev) <- gsub(colnames(dt_wide_prev), pattern = 'perc_nonzero\\.', replacement = '')

      # Set thresholds as a percentage of the maximum values in PP and GP
      pp_threshold_prev <- 0.8 * max(dt_wide_prev$PP)  # 80% of the max PP
      gp_threshold_prev <- 0.8 * max(dt_wide_prev$GP)  # 80% of the max GP

      # Filter rows where PP and GP exceed the thresholds
      filtered_df_prev <- dt_wide_prev[dt_wide_prev$PP > pp_threshold_prev | dt_wide_prev$GP > gp_threshold_prev, ]

      # Find the best row based on the criteria
      print(paste0('Prevalence: Optimal ',param_bnm, '_in_', tissue_area))
      if (nrow(filtered_df_prev) > 0) {
        best_row_prev <- filtered_df_prev[which.max(filtered_df_prev$PP + filtered_df_prev$GP - filtered_df_prev$IP), 'key']
        print(best_row_prev)
        cat('\n')
      } else {
        print("None meet the criteria.")
        cat('\n')
      }

      ########################################
      #### bar plots: distribution of  PP, IP, and GP phenotypes
      #### by individual detection rates
      ########################################

      ## ================
      ## load and combined cell density data from all input folders
      ## ================
      cell_all_cutoff <- c()
      for (ff in 1:length(folder_nms)){
        if(!file.exists((file.path(root_dir,folder_nms[ff],pheno_dir_path,pheno_count_fnm))))next
        setwd(file.path(root_dir,folder_nms[ff],pheno_dir_path))
        temp <- get(load(file=pheno_count_fnm))
        temp_rowids_area <- grep(x=temp$cell_ids, pattern = paste0('_',tissue_area,'$'))
        temp <- temp[temp_rowids_area,]
        temp$key <- folder_nms[ff]
        #cat(nrow(temp), '\n')
        cell_all_cutoff <- rbind(cell_all_cutoff, temp)
      }

      ## ================
      ## compute cell counts
      ## ================
      cell_count_by_phenotypes <- cell_all_cutoff  %>% ungroup() %>%
        group_by(key) %>% count(phenotype, name = 'no_of_cells')

      ## ================
      ## stacked bar-plot
      ## ================
      cell_count_by_phenotypes$detection_rate <- gsub(x=cell_count_by_phenotypes$key, pattern = param_bnm, replacement = '')
      cell_count_by_phenotypes$phenotype <- factor(cell_count_by_phenotypes$phenotype,
                                                   levels = c('PP','IP','GP'))
      pl2 <- ggplot(cell_count_by_phenotypes, aes(fill=phenotype, y=no_of_cells, x=detection_rate)) +
        geom_bar(position="stack", stat="identity")+
        theme_bw()+
        theme(axis.text = element_text(color='black',size = 11),
              axis.title = element_text(color='black',size=11),
              panel.border = element_rect(color='black',size=0.5),
              panel.grid = element_blank(),
              legend.spacing.x = unit(0.2, 'cm'),
              legend.title =  element_text(size = 11),
              legend.text = element_text(size = 11),
              legend.position = 'bottom',
              legend.direction = 'horizontal')+
        ylab('# cells')+ xlab(param_bnm)+
        labs(fill = "no. of cells")+
        guides(fill = guide_legend(title = "phenotypes", title.position = "left"))
      setwd(output_dir)
      ggsave(plot = pl2, width = plot_width, height = plot_height,
             filename = paste0(paste0('cellcounts_by',param_bnm, '_in_', tissue_area,'.pdf')))

      ## ==============
      # Compute %cell grouped by 'key'
      ## ==============
      cell_count_by_phenotypes <- cell_count_by_phenotypes %>%
        group_by(key) %>%
        mutate(percent_cell = (no_of_cells / sum(no_of_cells)) * 100)
      save(cell_count_by_phenotypes,
           file = file.path(output_dir,paste0('cellcounts_by',param_bnm, '_in_', tissue_area,'.RData')))

      ## ==============
      ## Identify the optimal setting
      ## ==============
      cell_count_by_phenotypes$no_of_cells <- NULL
      cell_count_by_phenotypes$detection_rate  <- NULL
      cell_count_by_phenotypes <- as.data.frame(cell_count_by_phenotypes)
      dt_wide <- reshape(cell_count_by_phenotypes, idvar = c("key"),
                         timevar = "phenotype", direction = "wide")
      ## replace NAs with Zeros
      dt_wide[, names(dt_wide) != "key"] <- lapply(dt_wide[, names(dt_wide) != "key"], function(x) ifelse(is.na(x), 0, x))
      ## clean column names
      colnames(dt_wide) <- gsub(colnames(dt_wide), pattern = 'percent_cell\\.', replacement = '')

      # Set thresholds as a percentage of the maximum values in PP and GP
      pp_threshold <- 0.8 * max(dt_wide$PP)  # 80% of the max PP
      gp_threshold <- 0.8 * max(dt_wide$GP)  # 80% of the max GP

      # Filter rows where PP and GP exceed the thresholds
      filtered_df <- dt_wide[dt_wide$PP > pp_threshold | dt_wide$GP > gp_threshold, ]

      # Find the best row based on the criteria
      print(paste0('%Cell: Optimal ',param_bnm, '_in_', tissue_area))
      if (nrow(filtered_df) > 0) {

        if(pp_threshold< 0.2*gp_threshold){
          best_row<-filtered_df[which.max(filtered_df$PP),'key']
        }else if(gp_threshold< 0.2*pp_threshold){
          best_row<-filtered_df[which.max(filtered_df$GP),'key']
        }else{
          best_row <- filtered_df[which.max(filtered_df$PP + filtered_df$GP - filtered_df$IP), 'key']
        }
        print(best_row)
        cat('\n')
      } else {
        print("None meet the criteria.")
        cat('\n')
      }
    }# end all tissue_area
  }

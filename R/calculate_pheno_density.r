#' Compute cell phenotype density
#'
#' Compute cell phenotype density
#' @param area_data A data.frame containing the area data (e.g. mm2) of each core (in rows, core_ids are
#' in row.names); the tissue categories (i.e. intraepithelial, stroma, overall) are represented in columns.
#' Note that matching annotation of tissue categories should be used in \code{area_data} and \code{phenotype_data}
#' identified in Step 3 previously. Note that "overall" naming is hard coded.
#' @param pheno_fnm A string indicating .RData object filename of the PIPA phenotype assignment
#' @param input_dir A file.path string of input directory storing the PIPA phenotype assignment of all cells of
#' interest
#' @param output_dir A file.path string of output directory
#' @import dplyr
#' @importFrom tidyr gather
#' @export
calculate_pheno_density <- function(area_data =NULL, input_dir=NULL, output_dir=NULL,
                              pheno_fnm='cell_phenotype.RData'){

  tumor_ids <- core_ids <- phenotype <- NULL

  ## ======================
  ## load cell phenotype data
  ## ======================
  setwd(input_dir)
  phenotype_data <- get(load(pheno_fnm))
  area <- 'overall'
  phenotype_data$area <- area

  ## ======================
  ## split by area
  ## ======================
  cellids_temp <- sapply(strsplit(x=phenotype_data$cell_ids, split = 'XY') , '[[', 2)
  possibleError <- tryCatch(
    area <- sapply(strsplit(x=cellids_temp, split = '_') , '[[', 3),
    error=function(e) e
  )
  if(!inherits(possibleError, "error")){
    phenotype_data$area <- area
  }
  cat('#cells in areas: ',  names(table(phenotype_data$area)), '\n')
  cat('#cells in areas: ',  table(phenotype_data$area), '\n')

  ## ======================
  ## phenotype count by area
  ## ======================
  area_core_counts <- phenotype_data %>% group_by(phenotype, area, core_ids,tumor_ids) %>% count(phenotype)

  ## ======================
  ## overall phenotype count
  ## ======================
  overall_core_counts <- phenotype_data %>% group_by(phenotype, core_ids,tumor_ids) %>% count(phenotype) %>%
    mutate(area = 'overall')
  ## ======================
  ## combine cell counts from all tissue cat/areas
  ## ======================
  core_counts <- rbind(area_core_counts, overall_core_counts)

  ## ======================
  ## area_data check: types of area/tissue segmentation
  ## ======================
  unique_area_in_cellcount <- unique(core_counts$area)
  missing_area <- unique_area_in_cellcount[!unique_area_in_cellcount %in% colnames(area_data)]
  if(length(missing_area) > 0) {
    warning(paste(missing_area, collapse = ','), ' area data is missing!\n')
    core_counts <- core_counts[core_counts$area != missing_area,]
  }

  ## ======================
  ## area_data check: remove area data not in phenotype data
  ## ======================
  area_data <- area_data[, colnames(area_data) %in% unique_area_in_cellcount, drop=FALSE]

  ## ======================
  ## area_data check: missing tumor_ids
  ## ======================
  uCore <- unique(core_counts$core_ids )
  missing_coreids <- uCore[! uCore %in% rownames(area_data)]
  if(length(missing_coreids) > 0){
    warning(': no area data provided!\n',paste0(missing_coreids, collapse = ','))
  }

  ## ======================
  ## convert area_data to long format: account for cores with 0 cell count but nonzero area which
  # should be considered as 0 density instead of missing value
  ## ======================
  pheno_by_area <- unique(core_counts[,c("phenotype","area")])
  tumor_core_ids_map <- unique(core_counts[,c("tumor_ids","core_ids")])

  area_data$core_ids <- rownames(area_data)
  area_colids <- setdiff(1:ncol(area_data), grep(x=colnames(area_data), pattern = 'core_ids'))
  area_data_long <- area_data %>% gather(area, area_value, area_colids)

  area_data_long <- merge(area_data_long, pheno_by_area, by=c("area"), all=TRUE)
  area_data_long <- merge(area_data_long, tumor_core_ids_map, by='core_ids')
  ## ======================
  ## merge area_data with phenotype count by core_ids
  ## ======================
  core_counts <- merge(core_counts, area_data_long, by = c('area','core_ids','phenotype','tumor_ids'),
                       all = TRUE)

  ## ======================
  ## area_data check: area == 0 but count >0
  ## ======================
  test <- core_counts %>% filter(!(is.na(n)) & area_value==0)
  if(nrow(test) > 0){
    setwd(output_dir)
    write.csv(test, file = 'error_zero_area.csv', row.names = FALSE)
    stop('area==0 for nonzero counts! check output file: error_zero_area.csv \n')
  }
  ## ======================
  ## exclude: area == 0 but count == 0
  ## ======================
  row_ids <- which( (is.na(core_counts$n)) & (core_counts$area_value==0) )
  if(length(row_ids)>0)core_counts <- core_counts[-row_ids, ]

  ## ======================
  ## exclude: area == NA
  ## ======================
  if(sum(is.na(core_counts$area_value)) > 0) {
    core_rm <- core_counts[(is.na(core_counts$area_value)),'core_ids']
    warning(paste0(core_rm, collapse = ','), '\nremoved!\n')
  }
  core_counts <- core_counts[!(is.na(core_counts$area_value)),]
  ## ======================
  ## replacing NA count with 0
  ## ======================
  core_counts$n[is.na(core_counts$n)] <- 0

  ## ======================
  ## calculate density: divide phenotype count by area
  ## ======================
  density <- area_value <- NULL
  core_counts <- core_counts %>% mutate(density= n/area_value)

  ## ======================
  ## calculate tumor-average cell density
  ## ======================
  tumor_avg_density <- core_counts %>% group_by(phenotype, area, tumor_ids) %>% summarise(density=mean(density, na.rm=TRUE))

  ## ======================
  ## saving
  ## ======================
  output_fnm <- gsub(x=pheno_fnm, pattern = '.RData', replacement = '_density.RData')
  setwd(output_dir)
  save(tumor_avg_density, file=output_fnm)


}#end function

#' Cell phenotype assignment
#'
#' Cell phenotype assignment for all cell of interest, using kNN
#'
#' @param cell_data A data.frame contains cells in rows and features in columns; cell_ids in row.names
#' the column 'tumor_ids' are the tumor identifiers; the column 'core_ids' are the core identifiers;
#' the column 'cell_ids' are the cell identifiers containing core_ids+coordinates
#' @param log10_featureNms A vector of characters of feature names whose values to be log10-transformed prior to scaling (z-score)
#' @param knn_neighbors An integers of the number of neighbors to be used \code{kNN}
#' @param selected_feature_dir A file.path string of input directory for selected features
#' (output from {meta_cluster})
#' @param seed An integer of seed for reproducibility
#' @param model_fnm A character string indicating the file name for the R object results output from
#' \code{meta_cluster}
#' @param min_fraction A fraction (0 to 1) indicating the minimum abundance of phenotype (among all cells of
#' interest) to be included in the differential analysis (default: 0.01 i.e. 1 percent)
#' @param cell_model_input_dir A file.path string of output directory that the cell model (i.e. meta cluster assignment)
#' to be loaded from
#' @param output_dir A file.path string of output directory that the resultant all cell phenotype assignment is saved
#' @param save_intermediate A boolean indicates the saving option for cell assignment of individual cell batches
#' @export
#' @importFrom DMwR kNN
#' @import dplyr
#'
cell_phenotype_assignment <- function(cell_data=NULL,log10_featureNms=NULL,
                                      selected_feature_dir=NULL,
                                      knn_neighbors= 100, seed=999,
                                      model_fnm="final_meta_clusters_by_prog.RData",
                                      cell_model_input_dir=NULL, output_dir=NULL,
                                      save_intermediate=FALSE,
                                      min_fraction=0.001){
  ## ======================
  ## parameters documentation
  ## ======================
  out_fnm <- file.path(output_dir,'cell_phenotype_assignment_log.txt')
  capture.output(file = out_fnm,append = FALSE,
                 cat('\n\n',format(Sys.time(), "%a %b %d %X %Y"),
                     '\nknn_neighbors=',knn_neighbors,
                     '\nmin_fraction=',min_fraction,
                     '\nmodel_fnm=',model_fnm,
                     '\nlog10_featureNms=' , log10_featureNms))
  ## ======================
  ## construct cell/core/tumor identifier data.frame:
  ## ======================
  ID <- cell_data[,c("cell_ids","core_ids","tumor_ids")]
  #ID <- unique(ID)
  rownames(cell_data) <- cell_data$cell_ids
  cell_data[,c("cell_ids","core_ids","tumor_ids")] <- NULL

  ## ======================
  ## load cell phenotyping model
  ## ======================
  setwd(file.path(cell_model_input_dir))
  meta_clusters <- get(load(model_fnm))

  ## ======================
  ## load user selected (cut-off) prognostic features
  ## ======================
  setwd(selected_feature_dir)
  selected_features <- read.csv(file = file.path('selectedCellFeatures.csv'), row.names = NULL,
                                as.is = TRUE)

  ## ======================
  ## clean up the (aggregate summary) feature names
  ## ======================
  selected_features <- selected_features[,1]
  capture.output(file = out_fnm,
                 cat('selected features:: ', paste0(selected_features, collapse = ', '), '\n'),
                 append = TRUE)

  cat('selected features:: ', paste0(selected_features, collapse = ', '), '\n')

  ## ======================
  ## filter cell data by selected features
  ## ======================
  cell_data <- cell_data[, colnames(cell_data) %in% selected_features, drop=FALSE]
  stopifnot(ncol(cell_data)==length(selected_features))

  ## ======================
  ## check if single cell feature is remaind
  ## ======================
  if(ncol(cell_data)==1){
    warning('ONLY 1 cell feature remains!')
    return(-1)
  }

  ## ======================
  ## identify training vs testing cells
  ## ======================
  ID$category <- ifelse(ID$cell_ids %in% meta_clusters$cell_ids, 'training', 'testing')
  capture.output(file = out_fnm,
                 cat('#cells in cell phenotyping model=', sum(ID$category=='training'), '\n'),
                 append = TRUE)
  capture.output(file = out_fnm,
                 cat('#cells for phenotype assignment=', sum(ID$category=='testing'), '\n'),
                 append = TRUE)

  cat('#cells in cell phenotyping model=', sum(ID$category=='training'), '\n')
  cat('#cells for phenotype assignment=', sum(ID$category=='testing'), '\n')

  ## ======================
  ## cleanup missing values
  ## ======================
  capture.output(file = out_fnm,
                 cat('#cells: before removing NAs=',nrow(cell_data), '\n'),
                 append = TRUE)
  cat('#cells: before removing NAs=',nrow(cell_data), '\n')
  cell_data <- cell_data[complete.cases(cell_data),, drop=FALSE]
  cell_data <- as.data.frame(cell_data)

  capture.output(file = out_fnm,
                 cat('#cells: after removing NAs=',nrow(cell_data), '\n'),
                 append = TRUE)
  cat('#cells: after removing NAs=',nrow(cell_data), '\n')

  ## ======================
  ## log10 transformation
  ## ======================
  for (cc in 1:ncol(cell_data)){
    feature_cc <- colnames(cell_data)[cc]
    if(feature_cc %in% log10_featureNms){
      cat('LOG10:', feature_cc, '\n')
      cell_data[,cc] <-log10_transform(cell_data[,cc], method = 'log10')
    }
  }

  ## ======================
  ## extract training cell data
  ## ======================
  cell_data_train <- cell_data[rownames(cell_data) %in% ID$cell_ids[ID$category=='training'], ,drop=FALSE]
  ## ======================
  ## calculate z-score of training cell data
  ## ======================
  cell_data_train <- scale(cell_data_train)

  ## ======================
  ## checking for invariant features in training cell data
  ## ======================
  invariant_features <- which(attr(cell_data_train,"scaled:scale")==0)
  invariant_features <- names(invariant_features)
  if(length(invariant_features)>0) warning('invariant features detected: ', invariant_features, '\n')

  ## ======================
  ## calculate z-score of all other cells based on the training model
  ## ======================
  cell_data_test <- cell_data[rownames(cell_data) %in% ID$cell_ids[ID$category=='testing'],, drop=FALSE]

  cell_data_test <- scale(cell_data_test, center = attr(cell_data_train,"scaled:center"),
                          scale = attr(cell_data_train, "scaled:scale"))

  ## ======================
  ## cleanup missing values
  ## ======================
  capture.output(file = out_fnm,
                 cat('#testing cells: before removing NAs=',nrow(cell_data_test), '\n'),
                 append = TRUE)
  cat('#testing cells: before removing NAs=',nrow(cell_data_test), '\n')
  cell_data_test <- cell_data_test[complete.cases(cell_data_test),,drop=FALSE]
  cell_data_test <- as.data.frame(cell_data_test)
  capture.output(file = out_fnm,
                 cat('#testing cells: after removing NAs=',nrow(cell_data_test), '\n'),
                 append = TRUE)
  cat('#testing cells: after removing NAs=',nrow(cell_data_test), '\n')

  ## ======================
  ## prepare class labels
  ## ======================
  phenotype<-0
  cell_data_test <- cbind(cell_data_test, phenotype)
  cell_data_train <- merge(meta_clusters[,c("cell_ids","phenotype")], cell_data_train,
                           by.x='cell_ids',by.y='row.names')
  rownames(cell_data_train) <- cell_data_train$cell_ids
  cell_data_train$cell_ids <- NULL
  ## ======================
  ## kNN cell phenotype assignment
  ## ======================
  cat('cell phenotype names (training model) : ', names(table(meta_clusters$phenotype)), '\n')
  cat('cell phenotype size (training model) : ', table(meta_clusters$phenotype), '\n')
  cat('kNN neighnorsize : ', knn_neighbors, '\n')

  capture.output(file = out_fnm,
                 cat('cell phenotype names (training model) : ', names(table(meta_clusters$phenotype)), '\n'),
                 append = TRUE)
  capture.output(file = out_fnm,
                 cat('cell phenotype size (training model) : ', table(meta_clusters$phenotype), '\n'),
                 append = TRUE)
  capture.output(file = out_fnm,
                 cat('kNN neighnorsize : ', knn_neighbors, '\n'),
                 append = TRUE)
  ## ======================
  ## randomize ordering
  ## ======================
  set.seed(seed)
  cell_data_train <- cell_data_train[sample(x = 1:nrow(cell_data_train), size = nrow(cell_data_train)),]
  ## ======================
  ## testing different k-neighbors
  ## ======================
  capture.output(file = out_fnm, cat('\n'),append = TRUE)
  capture.output(file = out_fnm, cat('kNN: k-neighbors = ', knn_neighbors, '\n'),
                 append = TRUE)


  if(nrow(cell_data_test)>300000){
    no_batch <- ceiling(nrow(cell_data_test)/300000)

    cell_data_test_id <- data.frame(id=c(1:nrow(cell_data_test)),stringsAsFactors = FALSE)
    cell_data_test_id$cut <- as.numeric(cut_number(cell_data_test_id$id,no_batch))
    nn5 <- c()
    for(bb in 1:no_batch){
      cell_data_test_sub <- cell_data_test[cell_data_test_id$id[cell_data_test_id$cut==bb], ]

      set.seed(seed)
      possibleError <- tryCatch(
        nn5_sub <- kNN(phenotype ~ .,cell_data_train,cell_data_test_sub,norm=FALSE,
                       k=knn_neighbors, prob=TRUE),
        error=function(e) e
      )
      if(inherits(possibleError, "error")) next
      setwd(output_dir)
      if(save_intermediate){
        save(nn5_sub, cell_data_test_sub, file = paste0('batch',bb,'_kNN_assignment.RData'))
      }
      nn5 <- c(nn5, as.character(nn5_sub))
    }#end test data no_batch
  }else{
    set.seed(seed)
    possibleError <- tryCatch(
      nn5 <- kNN(phenotype ~ .,cell_data_train,cell_data_test,norm=FALSE,k=knn_neighbors, prob=TRUE),
      error=function(e) e
    )
    if(inherits(possibleError, "error")) next
  }#end if (nrow(cell_data_test)>300000)


  capture.output(file = out_fnm,
                 cat('other cells assignemt to : ', names(table(nn5)), '\n'), append = TRUE)
  capture.output(file = out_fnm,
                 cat('other cells assignemt freq : ', table(nn5), '\n'), append = TRUE)


  capture.output(file = out_fnm,
                 cat('knn_neighbors=',knn_neighbors,',
                                   prob=', sum(attr(nn5, 'prob')), '\n'), append = TRUE)
  #### combind training & testing cell data
  phenotype_data <- data.frame(cell_ids= rownames(cell_data_test), phenotype= as.character(nn5), stringsAsFactors = FALSE)
  phenotype_data$core_ids <- sapply(strsplit(x=phenotype_data$cell_ids, split = 'XY') , '[[', 1)
  phenotype_data <- rbind(phenotype_data, meta_clusters)

  #### append tumor_ids to results
  phenotype_data <- merge(phenotype_data, ID[,c("cell_ids","tumor_ids")], by='cell_ids')

  #### saving
  setwd(output_dir)
  save(phenotype_data, file= paste0('raw_cell_phenotype.RData'))

  #### counting cells in each phenotype
  freq <- as.matrix(table(phenotype_data$phenotype))
  freq <- data.frame(phenotype =rownames(freq), cellcounts=freq[,1])
  fractions <- prop.table(table(phenotype_data$phenotype))
  freq$fractions <- fractions
  write.csv(freq, file = paste0('cellcounts_per_pheno.csv'),
            row.names = FALSE)

  #### excluding rare cell phenotype
  rare_pheno <- names(fractions)[which(fractions < min_fraction)]
  phenotype_data <- phenotype_data[!phenotype_data$phenotype %in% rare_pheno, ]
  #### saving
  setwd(output_dir)
  save(phenotype_data, file= paste0('cell_phenotype.RData'))


  #### counting cells in each phenotype grouped by areas (tumor vs stromal)
  cellids_temp <- sapply(strsplit(x=phenotype_data$cell_ids, split = 'XY') , '[[', 2)
  possibleError <- tryCatch(
    area <- sapply(strsplit(x=cellids_temp, split = '_') , '[[', 3),
    error=function(e) e
  )
  if(inherits(possibleError, "error")) stop()

  phenotype_data$area <- area
  cellCounts_byArea <- phenotype_data %>%  group_by(area) %>% count(phenotype, name = 'cellcounts')
  fractions <- prop.table(cellCounts_byArea$cellcounts)
  cellCounts_byArea$fractions <- fractions
  setwd(output_dir)
  write.csv(cellCounts_byArea, file = paste0('cellcounts_byPheno_byArea.csv'),
            row.names = FALSE)

  return(0)
}# end function

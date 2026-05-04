#' Refine phenotypes through merging PP or GP phenotypes based on prognostic significance
#'
#' Refine phenotypes through merging PP or GP phenotypes based on prognostic significance
#' (Step 1) Keep phenotypes show Ptrend < uni_P_cutoff in univariable Cox PH analysis
#' (Step 2) Keep phenotypes show association (i.e., non-zero coefficient) in a multivariable Cox
#' PH model which includes phenotypes remains from Step 1; separate models for PP and GP phenotypes.
#'
#' @param surv_data A data.frame: tumors in rows (tumor_ids in row.names), columns must include: "time", "cens";
#' if there are additional columns, they are treated as confounders
#' @param output_dir A file.path string of output directory
#' @param KM_xlab A character string for labelling the x-axis in Kaplan-Meier curves
#' @param preMerge_density_fnm A character string of input file name for cell density data before merging
#' @param preMerge_density_dir A file.path string of input directory storing the phenotype cell density data
#' before merging
#' @param preMerge_pheno_fnm A character string of input file name for cell phenotype data before merging
#' @param preMerge_pheno_dir A file.path string of input directory storing the cell phenotype data
#' before merging
#' @param uni_Cox_fnm A character string of input file name for univariable Cox PH analysis data
#' @param uni_Cox_dir A file.path string of input directory storing the univariable Cox PH analysis data
#' @param uni_P_cutoff A numeric value indicating the p-value cutoff for selecting phenotypes
#' @import dplyr survival survminer glmnet
#' @export
refine_pheno_byProg <- function(surv_data= NULL,  output_dir=NULL,
                            KM_xlab='Survivals',
                            preMerge_density_dir=NULL,preMerge_density_fnm = 'cell_phenotype_density.RData',
                            preMerge_pheno_dir=NULL,preMerge_pheno_fnm='cell_phenotype.RData',
                            uni_Cox_dir=NULL, uni_Cox_fnm='uniCoxPH_summary.txt',
                            uni_P_cutoff=0.005){
  read.table <- P_val <- V1 <- area <- pheno_area <- density <- tumor_ids <- all_of <- cens <- time <- NULL
  ## ======================
  ## parameters documentation
  ## ======================
  setwd(output_dir)
  capture.output(file = 'merge_phenotype_log.txt',
                 cat('preMerge_density_fnm=',preMerge_density_fnm,
                     '\nuni_Cox_dir=',uni_Cox_dir,
                     '\buni_Cox_fnm=',uni_Cox_fnm,
                     '\nuni_P_cutoff=',uni_P_cutoff))
  ## ======================
  ## identify confounders from input surv data
  ## ======================
  covariates <- setdiff(colnames(surv_data), c('cens','time'))

  ## ======================
  ## filter phenotypes
  ## ======================
  setwd(uni_Cox_dir)
  uni_COX_res <- read.table(file = uni_Cox_fnm, sep = ':', as.is = TRUE)

  uni_COX_res$P_val <- sapply(strsplit(x=uni_COX_res$V2, split = '='),"[[",2)
  uni_COX_res$P_val <- sapply(strsplit(x=uni_COX_res$P_val, split = ','),"[[",1)
  uni_COX_res$P_val <- as.numeric(uni_COX_res$P_val)
  uni_COX_res$V1 <- gsub(x=uni_COX_res$V1, pattern = ' ', replacement = '')

  uni_COX_res <- uni_COX_res %>% filter(P_val < uni_P_cutoff) %>% select(V1)
  uni_COX_res <- uni_COX_res %>% filter(grepl(V1, pattern = '_T|_S'))
  if(nrow(uni_COX_res)==0){
    warning("None of the phenotypes Ptrend is significant!")
    return(-1)
  }
  ## ======================
  ## load density data
  ## ======================
  setwd(preMerge_density_dir)
  tumor_avg_density <- get(load(file = preMerge_density_fnm))

  ## ======================
  ## cleanup density data
  ## ======================
  ## exclude other & overall area
  tumor_avg_density <- tumor_avg_density %>% filter(!area %in% c('overall', 'other'))

  ## rename area
  tumor_avg_density$pheno_area <- paste0(tumor_avg_density$phenotype, '_', tumor_avg_density$area)
  tumor_avg_density[, c("area","phenotype")] <- NULL
  tumor_avg_density_wide <- tidyr::spread(tumor_avg_density,
                                          pheno_area, density)
  ## ======================
  ## filter density data by significant phenotypes
  ## ======================
  tumor_avg_density_wide <- tumor_avg_density_wide %>% select(tumor_ids,contains('_T'), contains('_S'))
  tumor_avg_density_wide <- tumor_avg_density_wide[, colnames(tumor_avg_density_wide) %in%
                                                     c('tumor_ids', uni_COX_res[,1])]

  ## ======================
  ## calculate density quartile categories
  ## ======================
  for(r in 1:nrow(uni_COX_res)){
    pa <- uni_COX_res[r, 1]
    #cat(paste0(pa, collapse = ' in '), '\t')

    temp <- tumor_avg_density_wide %>% select(starts_with(pa)) %>% pull(!!pa)
    if( sum(is.na(temp))==nrow(tumor_avg_density_wide)){
      warning(cat(paste0(pa, collapse = ' in '), 'MISSING density data!'))
      next
    }
    q <- quantile(temp, na.rm = TRUE)
    #### check if too many zeros
    min_cut <- min(q)
    if(sum(q==min_cut)>1){
      cat0 <- temp[temp==min_cut]
      cat0 <- cat0[!is.na(cat0)]
      remaining <- temp[temp>min_cut]
      q3<- quantile(remaining, probs = seq(from = 0, 1, by = 0.33), na.rm = TRUE)
      q3[1] <- min(remaining, na.rm = TRUE)*0.99
      q3[length(q3)] <- max(remaining, na.rm = TRUE)
      q3_final <- c(min_cut,q3)
      if(sum(duplicated(q3_final))>0){
        density_data_sub$density_cat <- NULL
        cat('Too many zero\n')
        next
      }
      qq<- cut(temp, breaks=q3_final,
               include.lowest=TRUE)
      stopifnot(table(qq)[1]==length(cat0[!is.na(cat0)]))

    }else{
      qq<- cut(temp, breaks=quantile(temp, probs=seq(0,1, by=0.25), na.rm=TRUE),
               include.lowest=TRUE)
    }
    #### simplify/rename density Quartiles
    qq <- plyr::mapvalues(qq, from = levels(qq),
                          to = paste0(c('Q1','Q2','Q3','Q4')))
    #cat(table(qq),'\n')
    tumor_avg_density_wide[,pa] <- qq
  }## end all phenotypes



  ## ======================
  ## merge cell density + patient surv data
  ## ======================
  data <- merge(tumor_avg_density_wide, surv_data, by.x ='tumor_ids', by.y='row.names')

  ## ======================
  ## phenotype selection based on a multivariable Cox PH model, using Lasso method
  ## separate models for PP and Gp phenotypes
  ## ======================
  ## remove surv time equals to zero
  data <- data[data$time>0,]
  merging_phenotypes_id <- NULL
  for(aa in c('^PP','^GP')){
    aaa <- gsub(x=aa, pattern = '\\^', replacement = '')
    dens_colids <- grep(x=colnames(data), pattern = aa, value = TRUE)
    if(length(dens_colids)==0){
      setwd(output_dir)
      out_fnm <-  paste0(aaa,'.txt')
      capture.output(cat('no passing phenotype!'), file = out_fnm)
      next
    }

    ## set up ordinal variable
    data_ordinal <- data %>% select(all_of(dens_colids))
    data_ordinal <- data_ordinal %>% transmute_all(as.numeric)
    data_ordinal$tumor_ids <- data$tumor_ids

    ## append surv data
    data_ordinal <- data_ordinal %>% left_join(data %>% select(tumor_ids, all_of(covariates), cens, time), by='tumor_ids')

    ## build Cox PH model include only cell density categories
    x <- model.matrix(as.formula(paste( "~",paste0(dens_colids, collapse = "+"))), data_ordinal)
    y <- survival::Surv(data_ordinal$time, data_ordinal$cens)

    ## lamdba optimization by CV
    fit <- glmnet::glmnet(x, y, family="cox", alpha=1)
    set.seed(999)
    cv.fit <- glmnet::cv.glmnet(x, y, family="cox", alpha=1, nfolds=5)

    ## selected lambda
    min_lambda_id <- which(cv.fit$lambda==cv.fit$lambda.min)
    cv.fit$glmnet.fit$beta[1:nrow(cv.fit$glmnet.fit$beta), min_lambda_id]

    ## saving lasso coef
    setwd(output_dir)
    out_fnm <-  paste0(aaa,'_p',uni_P_cutoff,'.txt')
    capture.output(coef(cv.fit, s = "lambda.min"), file = out_fnm)

    remained_var <- names(which(abs(cv.fit$glmnet.fit$beta[,min_lambda_id])>1e-20))
    cat(remained_var, '\n')
    merging_phenotypes_id <- c(merging_phenotypes_id, remained_var)
  }#end phenotype progosis groups

  setwd(output_dir)
  save(merging_phenotypes_id, file = paste0('merging_phenotypes_id.RData'))


  ## ======================
  ## Post-PIPA:: re-assigning cells to merged phenotypes
  ## ======================
  ## load original phenotypes
  phenotype_data <- get(load(file= file.path(preMerge_pheno_dir,'cell_phenotype.RData')))
  merged_pheno <- phenotype_data$phenotype

  ## omit areas
  merging_phenotypes_id <- sapply(strsplit(x=merging_phenotypes_id, split = '_'), "[[", 1)

  ## merge PP cells
  PP_phenotypes <- grep(x=merging_phenotypes_id, pattern = '^PP', value = TRUE)
  merged_pheno[merged_pheno %in% PP_phenotypes] <- 'PP'

  ## merge GP cells
  GP_phenotypes <- grep(x=merging_phenotypes_id, pattern = '^GP', value = TRUE)
  merged_pheno[merged_pheno %in% GP_phenotypes] <- 'GP'

  ## merge remaining cells as IP cells
  merged_pheno[!merged_pheno %in% c('PP','GP')] <- 'IP'

  phenotype_data$phenotype <- merged_pheno
  cat("After merging:\n")
  cat((table(phenotype_data$phenotype)))

  ## saving
  setwd(output_dir)
  save(phenotype_data, file='cell_phenotype.RData')
  return(0)

}#end function

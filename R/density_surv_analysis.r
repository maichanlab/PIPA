#' Post-PIPA analysis - prognostic association analysis based on cell phenotype density
#'
#' Post-PIPA analysis - prognostic association analysis based on cell phenotype density
#' Patients/tumors are categorized into phenotype quartiles or zero+tertiles (if zero density cases >25% of total cases),
#' which are then used as (i) indicator variable in estimating Kaplan-Meier curves & log-rank test,
#' or (ii) as ordinal variable in Cox PH analysis
#' @param surv_data A data.frame: tumors in rows (tumor_ids in row.names), columns must include: "time", "cens";
#' if there are additional columns, they are treated as confounders
#' @param input_dir A file.path string of input directory storing the phenotype cell density data
#' @param output_dir A file.path string of output directory
#' @param KM_xlab A character string for labelling the x-axis in Kaplan-Meier curves
#' @param density_fnm A character string of input file name for cell density data
#' @import dplyr survival survminer
#' @export
density_surv_analysis <- function(surv_data= NULL, input_dir=NULL, output_dir=NULL,
                                  KM_xlab='Survivals', density_fnm = NULL){

  ## ======================
  ## identify cell phenotype result files
  ## ======================
  setwd(input_dir)
  ## get user-selected phenotype data
  if(is.null(density_fnm)){
    density_fnm <- paste0('cell_phenotype_density.RData')
  }


  ## ======================
  ## identify confounders
  ## ======================
  confounders <- setdiff(colnames(surv_data), c('cens','time'))

  ## ======================
  ## complete cases with survival
  ## ======================
  surv_data <- surv_data[complete.cases(surv_data[,c("cens","time")]),]


  ## ======================
  ## create output subdirectory
  ## ======================
  output_subdir <- file.path(output_dir, paste0('density_surv'))
  dir.create(output_subdir)

  ptrend_subdir <- file.path(output_subdir, 'ptrend')
  dir.create(ptrend_subdir)
  cat_subdir <- file.path(output_subdir, 'cat')
  dir.create(cat_subdir)
  km_subdir <- file.path(output_subdir, 'KM')
  dir.create(km_subdir)

  ## ======================
  ## load cell phenotype data
  ## ======================
  setwd(input_dir)
  density_data <- get(load(density_fnm))

  ## ======================
  ## subsetting dens data by complete cases with survival
  ## ======================
  density_data <- density_data[density_data$tumor_ids %in% rownames(surv_data), ]

  ## ======================
  ## create summary txt for all phenotype-area tests
  ## ======================
  uni_ptrend_summary_fnm <- file.path(ptrend_subdir,'uniCoxPH_summary.txt')
  capture.output(cat('\n'), file = uni_ptrend_summary_fnm)

  uni_cat_summary_fnm <- file.path(cat_subdir,'uniCoxPH_summary.txt')
  capture.output(cat('\n'), file = uni_cat_summary_fnm)

  if(length(confounders)>0){
    multi_ptrend_summary_fnm <- file.path(ptrend_subdir,'multiCoxPH_summary.txt')
    capture.output(cat('\n'), file = multi_ptrend_summary_fnm)
    multi_cat_summary_fnm <- file.path(cat_subdir,'multiCoxPH_summary.txt')
    capture.output(cat('\n'), file = multi_cat_summary_fnm)
  }
  ## ======================
  ## loop through each phenotype-area pair
  ## ======================
  phenotype <- area <- NULL
  pheno_area_pairs <- unique(density_data[,c("phenotype","area")])

  for(r in 1:nrow(pheno_area_pairs)){
    pa <- pheno_area_pairs[r, ]
    cat(paste0(pa, collapse = ' in '), '\t')

    density_data_sub <- merge(density_data, pa, by=c("phenotype","area"))

    #### categorize features selected: quartiles
    temp <- density_data_sub$density
    if( sum(is.na(temp))==nrow(density_data_sub)){
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
    cat(paste0(c('#tumors in : Q1','Q2','Q3','Q4')),'\n')
    cat(table(qq),'\n')
    density_data_sub$density_cat <- qq
    setwd(output_subdir)
    write.csv(x = density_data_sub,
              file = paste0('dens_quartiles_',pa$phenotype,'_',pa$area,'.csv'), row.names = FALSE)
    #### merge surv_data & cluster results
    m <- merge(density_data_sub, surv_data, by.x='tumor_ids', by.y = 'row.names')
    m <- m %>% mutate(phenotype_area = paste0(phenotype, '_',area))

    #### =============================================
    #### cast density cat to ordinal var
    #### =============================================
    m_ordinal <- m
    m_ordinal$density_cat <- as.numeric(m_ordinal$density_cat)
    key <- unique(m_ordinal$phenotype_area)
    m_ordinal[,key] <- m_ordinal$density_cat

    #### multivariate Cox PH
    if(length(confounders)>0){
      formula_mult <- as.formula(paste("Surv(time, cens)~",key,"+",paste0(confounders,collapse = '+')))
      res.cox <- survival::coxph(formula_mult, data =  m_ordinal)
      fileNm <- file.path(paste0('multiCoxPH_',paste0(pa,collapse = '_'),'.txt'))
      capture.output(summary(res.cox), file = file.path(ptrend_subdir,fileNm))


      s<-summary(res.cox)$coefficients
      CI<-summary(res.cox)$conf.int
      p.val <- round(s[1,"Pr(>|z|)"],digits = 5)
      capture.output(cat(rownames(s)[1],': p-value=',p.val[1], ', CI: ', CI[1,c(3,4)],'\n'),
                     file = multi_ptrend_summary_fnm,append = TRUE)
    }
    #### univariate Cox PH
    formula <- as.formula(paste("Surv(time, cens)~",key))
    res.cox <- survival::coxph(formula, data =  m_ordinal)
    fileNm <- file.path(paste0('uniCoxPH_',paste0(pa,collapse = '_'),'.txt'))
    capture.output(summary(res.cox), file = file.path(ptrend_subdir,fileNm))

    s<-summary(res.cox)$coefficients
    CI<-summary(res.cox)$conf.int
    p.val <- round(s[1,"Pr(>|z|)"],digits = 5)
    capture.output(cat(rownames(s)[1],': p-value=',p.val[1],
                       ', CI: ', CI[1,c(3,4)],'\n'),
                   file = uni_ptrend_summary_fnm,append = TRUE)

    #### =============================================
    #### categorical density
    #### =============================================
    key <- unique(m$phenotype_area)
    m[,key] <- m$density_cat

    ## output data container var
    multi_HR_CI <- list(); ll <- 1
    #### multivariate Cox PH
    if(length(confounders)>0){
      formula_mult <- as.formula(paste("Surv(time, cens)~",key,"+",paste0(confounders,collapse = '+')))
      res.cox <- survival::coxph(formula_mult, data =  m)
      fileNm <- file.path(paste0('multiCoxPH_',paste0(pa,collapse = '_'),'.txt'))
      capture.output(summary(res.cox), file = file.path(cat_subdir,fileNm))

      s<-summary(res.cox)$coefficients
      CI<-summary(res.cox)$conf.int
      for(rr in (grep(x=rownames(s), pattern = paste0(pa$phenotype,'_',pa$area))) ){
        p.val <- round(s[rr,"Pr(>|z|)"],digits = 5)
        capture.output(cat(rownames(s)[rr],'(spl size=',nrow(m),'), p-value=',p.val,
                           ', HR=', CI[rr,1],
                           ', lower95=', CI[rr,3], ',upper95=',CI[rr,4], '\n'),
                       file = file.path(multi_cat_summary_fnm), append = TRUE)
      }

      ## output saving
      multi_HR_CI[[ll]] <- summary(res.cox)$conf.int
      names(multi_HR_CI)[ll] <- levels(m$key)[1]
      ll <- ll + 1
      fileNm1 <- gsub(x=fileNm, pattern = '.txt', replacement = '_HR_CI.RData')
      save(multi_HR_CI, file=file.path(cat_subdir, fileNm1))
    }
    #### univariate Cox PH
    ## output data container var
    uni_HR_CI <- list(); ll <- 1

    formula <- as.formula(paste("Surv(time, cens)~",key))
    res.cox <- survival::coxph(formula, data =  m)
    fileNm <- file.path(paste0('uniCoxPH_',paste0(pa,collapse = '_'),'.txt'))
    capture.output(summary(res.cox), file = file.path(cat_subdir,fileNm))

    s<-summary(res.cox)$coefficients
    CI<-summary(res.cox)$conf.int
    for(rr in 1:nrow(s)){
      p.val <- round(s[rr,"Pr(>|z|)"],digits = 5)
      capture.output(cat(rownames(s)[rr],'(spl size=',nrow(m),'), p-value=',p.val,
                         ', HR=', CI[rr,1],
                         ', lower95=', CI[rr,3], ',upper95=',CI[rr,4], '\n'),
                     file = file.path(uni_cat_summary_fnm), append = TRUE)}

    ## output saving
    uni_HR_CI[[ll]] <- summary(res.cox)$conf.int
    names(uni_HR_CI)[ll] <- levels(m$key)[1]
    ll <- ll + 1
    fileNm1 <- gsub(x=fileNm, pattern = '.txt', replacement = '_HR_CI.RData')
    save(uni_HR_CI, file=file.path(cat_subdir, fileNm1))
    #### =============================================
    #### Kaplan-Meier curves
    #### =============================================
    fit <- survfit(formula, data =  m)
    d <- data.frame(time = fit$time,
                    n.risk = fit$n.risk,
                    n.event = fit$n.event,
                    n.censor = fit$n.censor,
                    surv = fit$surv,
                    upper = fit$upper,
                    lower = fit$lower
    )
    fit <- surv_fit(formula, data = m)

    pl<-ggsurvplot(fit,ylim = c(0, 1),
                   pval = TRUE,
                   pval.size=8,
                   fontsize = c(6),
                   font.title = c(16, "bold"),
                   font.subtitle = c(15, "bold"),
                   font.caption = c(16, "plain"),
                   font.x = c(18, "bold"),
                   font.y = c(18, "bold"),
                   font.tickslab = c(18),
                   font.legend = c(16),
                   risk.table.y.text = FALSE,
                   #trim.strata.names = TRUE,
                   xlab=KM_xlab,
                   risk.table = TRUE, # Add risk table
                   risk.table.col = "strata", # Change risk table color by groups
                   #linetype = "strata", # Change line type by groups
                   conf.int = FALSE,
                   surv.median.line = "hv", # Specify median survival
                   ggtheme = theme_bw()) # Change ggplot2 theme

    pl$plot <- pl$plot +
      theme(legend.text = element_text(size = 14, color = "black", face = "bold"))+
      guides(color = guide_legend(size=10,nrow=2,byrow=TRUE))

    plot_fileNm <- file.path(paste0('KM_',paste0(pa,collapse = '_'),'.pdf'))
    setwd(km_subdir)
    pdf(plot_fileNm, onefile = FALSE)
    print(pl)
    dev.off()


  }# end pheno_area_pairs


}#end function

# PIPA
Prognosis-informed Phenotype Assignment (PIPA): A Novel Method for Unsupervised Discovery of Cell Phenotypes with Prognostic Significance in the Tumor Microenvironment

## Initialization
````
# Package installation
# devtools::install_github(SIgN-CI/PIPA')
package_root_dir <- 'C:/Desktop'
install.packages(file.path(package_root_dir,'PIPA.tar.gz'), repos = NULL, type = 'source')
library(PIPA)

# Set paths
root_dir <- "C:/Desktop/";
results_folder <- 'PIPA_analysis'
output_dir <- file.path(root_dir, results_folder)
dir.create(output_dir)

# Set parameters
min_cluster_size <- 5 # Recommend min_cluster_size=20 for statistical significance in survival analysis
````

## Step 0: Prepare data for feature selection
````
# prepare cell feature data
cell_data_sub <- PIPA::cell_data
cell_data_sub[,c("celltype","core_ids","cell_ids")] <- NULL
# prepare surv data
surv_data <- patient_data[, c("cens","time")]
````


## Step 1: feature selection
````
step1_dir <- output_dir
log10_featureNms=c('feature6','feature7')
featureSelection_wrapper(data=cell_data_sub, surv_data=surv_data,
                         log10_featureNms=log10_featureNms,  
                         #no_of_runs=Lasso_run_no,
                         output_dir=step1_dir,
                         min_cluster_size= min_cluster_size)
````

## build phenotyping model at each identified feature subsets at various detection rates (steps 2-4 & post-analysis viz)
````
# get folder names for each candidate feature subsets
feature_subset_fnms <- list.files(path = step1_dir, pattern = "^detection_rate")

# loop over individual candidate feature subsets
for (ff in feature_subset_fnms){
  subtype_subdir <- file.path(step1_dir, ff)
  if(!dir.exists(subtype_subdir))next
  
# Step 2: identify meta clusters via consensus clustering and merging 
model_dir <- file.path(subtype_subdir, 'step2_metaCluster')
dir.create(model_dir)
cell_data <- PIPA::cell_data
  meta_cluster(cell_data=cell_data,log10_featureNms=log10_featureNms,
               input_dir=subtype_subdir,output_dir=model_dir,
               batch_size=500, 
               consensus_rep=2,
               no_of_random_merging=5,
               max_batch=5,
               cell_batch_clustering = TRUE,
               cell_batch_merging=TRUE)
  
# Step 3:  cell phenotype assignment of all cells of interest based on the meta clusters
pheno_assignment_dir <- file.path(subtype_subdir, 'step3_pheno_assignment')
dir.create(pheno_assignment_dir)
cell_phenotype_assignment(cell_data=cell_data,log10_featureNms=log10_featureNms,selected_feature_dir=subtype_subdir,cell_model_input_dir=model_dir,output_dir=pheno_assignment_dir)
  
# Step 4: refine phenotypes by merging cell clusters that share a significant and consistent prognostic direction
  
# (i) calculate phenotype density
  calculate_pheno_density(area_data=area_data,                            input_dir=pheno_assignment_dir, output_dir=pheno_assignment_dir)
  
# (ii) perform cell density-based survival analysis
final_pheno_dir <- file.path(subtype_subdir, 'step4_final_pheno')
dir.create(final_pheno_dir)
colnames(patient_data)
density_surv_analysis(surv_data= patient_data,                         input_dir=pheno_assignment_dir, output_dir=final_pheno_dir,                         density_fnm= 'cell_phenotype_density.RData',                        KM_xlab= 'survivals')
  
# (iii) cell clusters merging based on prognostic significance and direction
cutoff_pval <- 0.05 # this is only for demonstration, please use cutoff 0.005
merged_pheno_dir <- file.path(final_pheno_dir, paste0('merged_pheno_P',cutoff_pval))
dir.create(merged_pheno_dir)
out <- refine_pheno_byProg(surv_data= patient_data,                              KM_xlab='survivals', output_dir=merged_pheno_dir,                             preMerge_density_dir =pheno_assignment_dir,                              preMerge_pheno_dir=pheno_assignment_dir,
                             uni_Cox_dir=file.path(final_pheno_dir,'density_surv/ptrend'), 
uni_Cox_fnm='uniCoxPH_summary.txt',uni_P_cutoff=cutoff_pval)
if(out==(-1))next
  
# Step 5 Post-PIPA analysis: survival analysis using PIPA-derived PP/IP/GP phenotypes
# calculate cell density of PIPA-derived PP/IP/GP phenotypes
calculate_pheno_density(area_data=area_data,                            input_dir=merged_pheno_dir, output_dir=merged_pheno_dir)
  
# (i) phenotype density-based survival analysis  density_surv_analysis(surv_data= patient_data,                          input_dir=merged_pheno_dir, output_dir=merged_pheno_dir,                         density_fnm= 'cell_phenotype_density.RData',                        KM_xlab= 'Survivals')
  
# (ii) composition-based survival analysis   composition_surv_analysis(surv_data= patient_data,                              density_data = raw_dens_data,                            min_cluster_size = min_cluster_size,                             phenotype_dir=merged_pheno_dir, output_dir=merged_pheno_dir,                            neighbor_size_byFrac= seq(0.05,0.25,by=0.15),                            KM_xlab='Survivals',area_split=TRUE)
  
# Step 6 Post-PIPA analysis: feature visualization
# retrieve PIPA-selected cell features
selected_feature <- read.csv(file = file.path(subtype_subdir,'selectedCellFeatures.csv'), as.is = TRUE)
  
# (i) feature heatmap
# subsetting cell data for PIPA-selected cell features
cell_data_sub <- PIPA::cell_data[, c("cell_ids","tumor_ids",selected_feature$features)]
Viz_feature_heatmap_byPheno(cell_data= cell_data_sub,                              pheno_dir=merged_pheno_dir, output_dir=merged_pheno_dir,                              pheno_levels = c('PP','IP','GP'),column_scale=TRUE)
  
# (ii) feature boxplot
# subsetting cell data for PIPA-selected cell features
cell_data_sub <- PIPA::cell_data[, c("cell_ids",selected_feature$features)]
Viz_feature_boxplot_byPheno(cell_data= cell_data_sub,plot_bnm='PIPA_features',outlier_removal=TRUE, pheno_levels = c('PP','IP','GP'),                              pheno_dir=merged_pheno_dir, output_dir=merged_pheno_dir,                              ncol = 4, plot_height = 16, plot_width = 15)
  
# (iii) cell umap
# subsetting cell data for PIPA-selected cell features
cell_data_sub <- PIPA::cell_data[, c("cell_ids",selected_feature$features)]
Viz_feature_UMAP(cell_data= cell_data_sub, downspl_size=100,                   pheno_dir=merged_pheno_dir, output_dir=merged_pheno_dir,                   log10_featureNms=log10_featureNms,                   plot_height = 8, plot_width = 8)

}# end all candidate feature subsets (at different detection rates)
````



## Optimize parameter settings: min. detection rates for optimal feature selection in 'featureSelection_patientGrouping'
````
param_bnm='detection_rate'
output_dir2 <- file.path(step1_dir, paste0("Viz_all_",param_bnm)); dir.create(output_dir2)
cutoff_pval <- 0.05 # this is only for demonstration, please use cutoff 0.005
pheno_dir_path <- paste0('step4_final_pheno/merged_pheno_P',cutoff_pval)
# (i) identify the optimal min. detection rate that maximizes the abundance and prevalence of detected cell phenotypes
Optimize_Params(root_dir=step1_dir,
                              output_dir= output_dir2,
                              pheno_dir_path=pheno_dir_path,
                              param_bnm=param_bnm,
                              tissue_area_list=c('T','S'))

# check
list.files(file.path(step1_dir, 'Viz_all_detection_rate'), pattern = '.RData')
cellcounts_in_S <- get(load(file.path(step1_dir, 'Viz_all_detection_rate','cellcounts_bydetection_rate_in_S.RData')))
prevalent_in_T <- get(load(file.path(step1_dir, 'Viz_all_detection_rate','prevalence_bydetection_rate_in_T.RData')))

# (ii) generate forest plots from Cox proportional hazards regression results computed using the density_surv_analysis function across varying detection thresholds
Viz_foresPlot_denQ_CoxPH(root_dir=step1_dir, output_dir= output_dir2,                         pheno_dir_path=pheno_dir_path,                         param_bnm='detection_rate',                         Cox_result_folderNm='density_surv/cat')
````


## Optimize parameter settings: batch size in 'meta_cluster' function
````
#  Collect data for selected batch sizes
batch_size_list <- c(300, 400, 500, 600)
## select the optimal detection rate
subtype_subdir <- file.path(step1_dir, 'detection_rate0.3')
if(!dir.exists(subtype_subdir))next
model_dir <- file.path(subtype_subdir, 'step2_metaCluster')
dir.create(model_dir)

# Perform Steps 2–5 followed by post-analysis to derive the final phenotypes and their associated survival outcomes
for(bb in batch_size_list){
output_model_dir <- file.path(model_dir,paste0('batch',bb)) 
dir.create(output_model_dir)
# Step 2: identify meta clusters via consensus clustering and merging 
cell_data <- PIPA::cell_data
  meta_cluster(cell_data=cell_data,log10_featureNms=log10_featureNms,
               input_dir=subtype_subdir,output_dir=output_model_dir,
               batch_size=bb, 
               consensus_rep=2,
               no_of_random_merging=5,
               max_batch=5,
               cell_batch_clustering = TRUE,
               cell_batch_merging=TRUE)
  
  
# Step 3:  cell phenotype assignment of all cells of interest based on the meta clusters
pheno_assignment_dir <- file.path(output_model_dir, 'step3_pheno_assignment')
dir.create(pheno_assignment_dir)
  cell_phenotype_assignment(cell_data=cell_data,log10_featureNms=log10_featureNms,                            selected_feature_dir=subtype_subdir,                            cell_model_input_dir=output_model_dir,                            output_dir=pheno_assignment_dir)

# Step 4: refine phenotypes by merging cell clusters that share a significant and consistent prognostic direction
# (i) calculate phenotype density
calculate_pheno_density(area_data=area_data,                            input_dir=pheno_assignment_dir, output_dir=pheno_assignment_dir)
# (ii) perform cell density-based survival analysis
final_pheno_dir <- file.path(output_model_dir, 'step4_final_pheno')
dir.create(final_pheno_dir)
colnames(patient_data)
density_surv_analysis(surv_data= patient_data,  
                        input_dir=pheno_assignment_dir, output_dir=final_pheno_dir,                         density_fnm= 'cell_phenotype_density.RData',                        KM_xlab= 'survivals')
# (iii) cell clusters merging based on prognostic significance and direction
cutoff_pval <- 0.05 # this is only for demonstration, please use cutoff 0.005
merged_pheno_dir <- file.path(final_pheno_dir, paste0('merged_pheno_P',cutoff_pval))
dir.create(merged_pheno_dir)
out <- refine_pheno_byProg(surv_data= patient_data, 
                                KM_xlab='survivals', output_dir=merged_pheno_dir,                                preMerge_density_dir =pheno_assignment_dir, 
preMerge_pheno_dir=pheno_assignment_dir,                                uni_Cox_dir=file.path(final_pheno_dir,'density_surv/ptrend'),                                uni_Cox_fnm='uniCoxPH_summary.txt',uni_P_cutoff=cutoff_pval)
if(out==(-1))next
  
# Step 5 Post-PIPA analysis: survival analysis using PIPA-derived PP/IP/GP phenotypes
# calculate cell density of PIPA-derived PP/IP/GP phenotypes
calculate_pheno_density(area_data=area_data,  
                          input_dir=merged_pheno_dir, output_dir=merged_pheno_dir)
  
# (i) phenotype density-based survival analysis  
  density_surv_analysis(surv_data= patient_data,  
                        input_dir=merged_pheno_dir, output_dir=merged_pheno_dir,                         density_fnm= 'cell_phenotype_density.RData',                        KM_xlab= 'Survivals')
  
#(ii) composition-based survival analysis   
  composition_surv_analysis(surv_data= patient_data,                              density_data = raw_dens_data,                            min_cluster_size = min_cluster_size, 
phenotype_dir=merged_pheno_dir, output_dir=merged_pheno_dir,                            neighbor_size_byFrac= seq(0.05,0.25,by=0.15),
KM_xlab='Survivals',area_split=TRUE)
  
# Step 6 Post-PIPA analysis: feature visualization
# retrieve PIPA-selected cell features
selected_feature <- read.csv(file = file.path(subtype_subdir,'selectedCellFeatures.csv'), as.is = TRUE)
  
# (i) feature heatmap
cell_data_sub <- PIPA::cell_data[, c("cell_ids","tumor_ids",selected_feature$features)]
Viz_feature_heatmap_byPheno(cell_data= cell_data_sub,                              pheno_dir=merged_pheno_dir, output_dir=merged_pheno_dir,                              pheno_levels = c('PP','IP','GP'),                              column_scale=TRUE)
# (ii) feature boxplot
cell_data_sub <- PIPA::cell_data[, c("cell_ids",selected_feature$features)]
Viz_feature_boxplot_byPheno(cell_data= cell_data_sub,plot_bnm='PIPA_features',outlier_removal=TRUE, pheno_levels = c('PP','IP','GP'),                              pheno_dir=merged_pheno_dir, output_dir=merged_pheno_dir,                              ncol = 4, plot_height = 16, plot_width = 15)
  
# (iii) cell umap
cell_data_sub <- PIPA::cell_data[, c("cell_ids",selected_feature$features)]
Viz_feature_UMAP(cell_data= cell_data_sub, downspl_size=100,              pheno_dir=merged_pheno_dir, output_dir=merged_pheno_dir,                   log10_featureNms=log10_featureNms,                   plot_height = 8, plot_width = 8)
}#end all batch size

# Identify optimal batch size for meta_cluster
param_bnm='batch'
output_dir2 <- file.path(model_dir, paste0("Viz_all_",param_bnm)); dir.create(output_dir2)
cutoff_pval <- 0.05 # this is only for demonstration, please use cutoff 0.005
pheno_dir_path <- paste0('step4_final_pheno/merged_pheno_P',cutoff_pval)
Optimize_Params(root_dir=model_dir,                              output_dir= output_dir2,                              pheno_dir_path=pheno_dir_path,                              param_bnm=param_bnm,                              tissue_area_list=c('T','S'))
````


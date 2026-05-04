#' Cell feature data.
#'
#' A dataset containing 20000 cell of interest,
#' 15793 cells in stromal region and 4207 cells in tumor region
#'
#' @format A data frame with 10000 rows and 3 identity columns + 19 cell feature columns: \describe{
#'   \item{feature1}{cell feature 1, numeric values} \item{feature2}{cell feature 2, numeric values}
#'   \item{feature3}{cell feature 3, numeric values} \item{feature4}{cell feature 4, numeric values}
#'   \item{feature5}{cell feature 5, numeric values} \item{feature6}{cell feature 6, numeric values}
#'   \item{feature7}{cell feature 7, numeric values} \item{feature8}{cell feature 8, numeric values}
#'   \item{feature9}{cell feature 9, numeric values} \item{feature10}{cell feature 10, numeric values}
#'   \item{feature11}{cell feature 11, numeric values} \item{feature12}{cell feature 12, numeric values}
#'   \item{feature13}{cell feature 13, numeric values} \item{feature14}{cell feature 14, numeric values}
#'   \item{feature15}{cell feature 15, numeric values} \item{feature16}{cell feature 16, numeric values}
#'   \item{feature17}{cell feature 17, numeric values} \item{feature18}{cell feature 18, numeric values}
#'   \item{feature19}{cell feature 19, numeric values}
#'   \item{cell_ids}{cell identifiers with XY locations, character strings}
#'   \item{tumor_ids}{tumor or patient identifiers, character strings}
#'   \item{core_ids}{core identifiers, character strings} }
"cell_data"
#' tissue core area data.
#'
#' A data.frame containing containing 200 rows for 100 cores comprising 100 cases/tumors
#'
#' @format A data.frame with 200 rows (core_ids in row.names) and 3 columns:
#' \describe{
#'   \item{S}{stromal area in mm2, positive numeric values}
#'   \item{T}{tumor intraepitheliao area in mm2, positive numeric values}
#'   \item{overall}{overall area in mm2, positive numeric values}
#'   }
"area_data"
#' Tumor denisty data.
#'
#' A data.frame containing containing 200 rows, 100 each in stromal and
#' tumor area (i.e. 100 tumors)
#'
#' @format A data.frame with 200 rows and 4 columns:
#' \describe{
#'   \item{tumor_ids}{tumor or patient identifiers, character strings}
#'   \item{area}{S for stromal and T for tumor intraepithelial region, character strings}
#'   \item{phenotype}{cell phenotype name, character strings}
#'   \item{density}{density values, positive numeric values}
#'   }
"raw_dens_data"
#' Clinical outcome and clinicopathological data.
#'
#' A data frame containing the survival and clinicopathological data for 100 tumors.
#'
#' @format A data frame containing 100 rows (i.e. cases/tumors) and 6 columns:
#'   \describe{ \item{time}{survival time, numeric} \item{cens}{censoring data,
#'   binary} \item{cov1}{covariate var1, factor}
#'   \item{cov2}{covariate var, factor} \item{cov3}{covariate var, factor}
#'   \item{cov4}{covariate var, factor}}
"patient_data"

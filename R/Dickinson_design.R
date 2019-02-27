#' @name Dickinson_design
#' @title Raw county-level variables for study 1 in Dickinson et al (2015)
#' @description
#'   Two approaches (interventions) are compared for increasing the "up-to-date" immunization rate
#'  in 19- to 35-month-old children. 16 counties in Colorado 1:1 are randomized to either a
#'  population-based approach or a practice-based approach. Ahead of randomization, several
#'  county-level variables are collected, and a subset of them are used for covariate constrained
#'  randomization. The continuous variable of average income is categorized to illustrate the
#'  use of cvcrand on multi-category variables. And the percentage in CIIS variable is truncated at 100%.
#' @docType data
#' @format A data frame with 16 rows and 7 variables:
#' \describe{
#'   \item{county}{the identification for the county}
#'   \item{location}{urban or rural}
#'   \item{inciis}{percentage of children aged 19-35 months in the Colorado Immunization Information System (CIIS)}
#'   \item{numberofchildrenages1935months}{number of children aged 19-35 months}
#'   \item{uptodateonimmunizations}{percentage of children already up-to-date on their immunization}
#'   \item{africanamerican}{percentage of African American}
#'   \item{hispanic}{percentage of Hispanic ethnicity}
#'   \item{income}{average income}
#'   \item{incomecat}{average income categorized into tertiles}
#'   \item{pediatricpracticetofamilymedicin}{pediatric practice-to-family medicine practice ratio}
#'   \item{communityhealthcenters}{number of community health centers}
#' }
#' @source \url{https://www.jabfm.org/content/28/5/663/tab-figures-data}
#' @references
#'   Dickinson, L. M., B. Beaty, C. Fox, W. Pace, W. P. Dickinson, C. Emsermann,
#'   and A. Kempe (2015): Pragmatic cluster randomized trials using covariate
#'   constrained randomization: A method for practice-based research networks (PBRNs).
#'   The Journal of the American Board of Family Medicine 28(5): 663-672
NULL

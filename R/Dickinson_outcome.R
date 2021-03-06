#' @name Dickinson_outcome
#' @title Simulated individual-level binary outcome and baseline variables for study 1 in Dickinson et al (2015)
#' @description
#'   At the end of the study, the researchers will have ascertained the outcome in the
#'   16 clusters. Suppose that the researchers were able to assess 300 children in each cluster.
#'   We simulated correlated outcome data at the individual level using a generalized linear
#'   mixed model (GLMM) to induce correlation by including a random effect. The intracluster
#'   correlation (ICC) was set to be 0.01, using the latent response definition provided in
#'   Eldrige et al. (2009). This is a reasonable value of the ICC for population health studies
#'   (Hannan et al. 1994). We simulated one data set, with the outcome data dependent on the
#'   county-level covariates used in the constrained randomization design and a positive treatment
#'   effect so that the practice-based intervention increases up-to-date immunization rates more
#'   than the community-based intervention. For each individual child, the outcome is equal to 1
#'   if he or she is up-to-date on immunizations and 0 otherwise.
#'
#'   Note that we still categorize the continuous variable of average income to illustrate
#'   the use of cvcrand on multi-category variables, and we truncated the percentage
#'   in CIIS variable at 100%.
#' @docType data
#' @format A data frame with 4800 rows and 7 variables:
#' \describe{
#'   \item{county}{the identification for the county}
#'   \item{location}{urban or rural}
#'   \item{inciis}{percentage of children ages 19-35 months in the Colorado Immunization Information System (CIIS)}
#'   \item{uptodateonimmunizations}{percentage of children already up-to-date on their immunization}
#'   \item{hispanic}{percentage of population that is Hispanic}
#'   \item{incomecat}{average income categorized into tertiles}
#'   \item{outcome}{the status of being up-to-date on immunizations}
#' }
#' @references
#'   Dickinson, L. M., B. Beaty, C. Fox, W. Pace, W. P. Dickinson, C. Emsermann,
#'   and A. Kempe (2015): Pragmatic cluster randomized trials using covariate
#'   constrained randomization: A method for practice-based research networks (PBRNs).
#'   The Journal of the American Board of Family Medicine 28(5): 663-672
#'
#'   Eldridge, S. M., Ukoumunne, O. C., & Carlin, J. B. (2009). The Intra Cluster
#'   Correlation Coefficient in Cluster Randomized Trials: A Review of Definitions.
#'   International Statistical Review, 77(3), 378-394.
#'
#'   Hannan, P. J., Murray, D. M., Jacobs Jr, D. R., & McGovern, P. G. (1994).
#'   Parameters to aid in the design and analysis of community trials: intraclass
#'   correlations from the Minnesota Heart Health Program. Epidemiology, 88-95. ISO 690
NULL

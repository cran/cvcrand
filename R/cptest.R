#' Clustered permutation test for cluster randomized trials
#' @param outcome a vector specifying the individual-level outcome.
#' @param clustername a vector specifying the identification variable of the cluster.
#' @param z  a data frame of covariates to be adjusted for in the permutation analysis.
#' @param cspacedatname  gives the path of the csv dataset containing the saved randomization space. This dataset contains the permutation matrix, as well as an indicator variable in the first column indicating which row of the permutation matrix was selected as the final scheme to be implemented in practice.
#' @param outcometype the type of regression model that should be run. Options are \code{"continuous"} for linear regression and \code{"binary"} for logistic regression.
#' @param categorical a vector specifying categorical (including binary) variables. This can be names of the columns or number indexes of columns, but cannot be both. Suppose there are \code{p} categories for a categorical variable, \code{cptest} function creates \code{(p-1)} dummy variables and drops the reference level if the variable is specified as a factor. Otherwise, the first level in the alphanumerical order will be dropped. If the user wants to specify a different level to drop for a \code{p}-level categorical variable, the user can create \code{p-1} dummy variables and these can instead be supplied as covariates to the \code{cptest} function. Then, the user needs to specify the dummy variables created themselves to be \code{categorical} when running \code{cptest}. In addition, the user could also set the variable as a factor with the specific reference level. The user must ensure that the same level of the categorical variable is excluded as was excluded when running \code{cvrall}, by coding the variables the same way as in the design phase. This is the only optional argument of the \code{cptest} function. All others are required.
#' @keywords cluster-randomized-trails clustered-permutation-test
#' @author Hengshi Yu <hengshi@umich.edu>, Fan Li <fan.f.li@yale.edu>, John A. Gallis <john.gallis@duke.edu>, Elizabeth L. Turner <liz.turner@duke.edu>
#' @description cptest performs a clustered permutation test on the individual-level outcome data for cluster
#' randomized trials (CRTs). The type of the outcome can be specified by the user to be \code{"continuous"} or
#' \code{"binary"}.
#'
#' Linear regression (for outcome type \code{"continuous"}) or logistic regression (for outcome type \code{"binary"}) is applied to the outcome regressed on covariates specified. Cluster residual means are computed. Within the constrained space,
#' the contrast statistic between the treatment and control arms is created from the randomization schemes and the cluster residual means. The permutation test is then conducted by comparing the contrast statistic for the scheme actually utilized to all other schemes in the constrained space.
#' @references
#' Gail, M.H., Mark, S.D., Carroll, R.J., Green, S.B. and Pee, D., 1996. On design considerations and randomization based inference for community intervention trials. Statistics in medicine, 15(11), pp.1069-1092.
#'
#' Li, F., Lokhnygina, Y., Murray, D.M., Heagerty, P.J. and DeLong, E.R., 2016. An evaluation of constrained randomization for the design and analysis of group randomized trials. Statistics in medicine, 35(10), pp.1565-1579.
#'
#' Li, F., Turner, E. L., Heagerty, P. J., Murray, D. M., Vollmer, W. M., & DeLong, E. R. (2017). An evaluation of constrained randomization for the design and analysis of group randomized trials with binary outcomes. Statistics in medicine, 36(24), 3791-3806.
#'
#' Gallis, J. A., Li, F., Yu, H., Turner, E. L. (In Press).  cvcrand and cptest: Efficient design and analysis of cluster randomized trials.  Stata Journal.
#'
#' Gallis, J. A., Li, Fl. Yu, H., Turner, E. L. (2017).  cvcrand and cptest: Efficient design and analysis of cluster randomized trials.  Stata Conference.  https://www.stata.com/meeting/baltimore17/slides/Baltimore17_Gallis.pdf.
#'
#' Dickinson, L. M., Beaty, B., Fox, C., Pace, W., Dickinson, W. P., Emsermann, C., & Kempe, A. (2015). Pragmatic cluster randomized trials using covariate constrained randomization: A method for practice-based research networks (PBRNs). The Journal of the American Board of Family Medicine, 28(5), 663-672.
#' @export
#' @examples
#'\dontrun{
#' Analysis_result <- cptest(outcome = Dickinson_outcome$outcome,
#'                           clustername = Dickinson_outcome$county,
#'                           z = data.frame(Dickinson_outcome[ , c("location", "inciis",
#'                               "uptodateonimmunizations", "hispanic", "incomecat")]),
#'                           cspacedatname = "dickinson_constrained.csv",
#'                           outcometype = "binary",
#'                           categorical = c("location","incomecat"))
#'}
#'
#' @return \code{FinalScheme} the final scheme in the permutation matrix
#' @return \code{pvalue} the p-value of the intervention effect from the clustered permutation test
#' @return \code{pvalue_statement} the statement about the p-value of the intervention effect from the clustered permutation test






cptest <- function(outcome, clustername, z = NULL, cspacedatname, outcometype, categorical = NULL){
    x <- z
  if(!is.null(z)){
       if (is.null(categorical)) { # if there are no categorical variables specified

   if (sum(apply(x, 2, function(x) length(unique(x[!is.na(x)]))) <= 1) >= 1) {
     ## check the variables with only one unique value, i.e. no variation

     warning(cat("Warning: each of columns", c(1 : dim(x)[2])[apply(x, 2, function(x) length(unique(x[!is.na(x)]))) <= 1], "has only less than or equal to one unique value and no variation."))
   }




   if (!setequal(c(1 : dim(x)[2])[apply(x, 2, function(x) length(unique(x[!is.na(x)]))) <= 4], categorical)) {
     ## check the categorical variables have been correctly specified

     warning(cat("Warning: each of columns", c(1 : dim(x)[2])[apply(x, 2, function(x) length(unique(x))) <= 4], "has less than or equal to 4 unique observations. Please check all the categorical variables have been correctly specified."))


   }


    } else {  # If there are categorical variables specified

    p <- length(categorical)

     if (is.character(categorical)) {  ## columns names based categorical variables

      if (sum(apply(x, 2, function(x) length(unique(x[!is.na(x)]))) <= 1) >= 1) {
         ## check the variables with only one unique value, i.e. no variation

         warning(cat("Warning: each of", colnames(x)[apply(x, 2, function(x) length(unique(x[!is.na(x)]))) <= 1], "has only less than or equal to one unique value and no variation."))
       }


      if (!setequal(colnames(x)[apply(x, 2, function(x) length(unique(x[!is.na(x)]))) <= 4], categorical)) {
        ## check the categorical variables have been correctly specified


        warning(cat("Warning: each of", setdiff(colnames(x)[apply(x, 2, function(x) length(unique(x))) <= 4], categorical), "has less than or equal to 4 unique observations. Please check all the categorical variables have been correctly specified."))
      }



       for(i in 1:p){

             x <- data.frame(x, model.matrix(~factor(x[ , categorical[i]]))[ , -1])
                           # add the dummy variables into the x matrix or data frame
                     }


        x <- x[ , -which(colnames(x) %in% categorical)]
                           # remove the initial categorical variables from x

          } else if (is.numeric(categorical)) { ## columns indexes based categorical variables

            if(sum(apply(x, 2, function(x) length(unique(x[!is.na(x)]))) <= 1) >= 1){
             ## check the variables with only one unique value, i.e. no variation

           warning(cat("Warning: each of columns", c(1:dim(x)[2])[apply(x, 2, function(x) length(unique(x[!is.na(x)]))) <= 1], "has only less than or equal to one unique value and no variation."))

          }


         if (!setequal(c(1:dim(x)[2])[apply(x, 2, function(x) length(unique(x[!is.na(x)]))) <= 4], categorical)) {
                           ## check the categorical variables have been correctly specified

       warning(cat("Warning: each of columns",setdiff(c(1:dim(x)[2])[apply(x, 2, function(x) length(unique(x))) <= 4], categorical), "has less than or equal to 4 unique observations, Please check all the categorical variables have been correctly specified."))


     }



         for (i in 1:p) {

               x <- data.frame(x, model.matrix(~factor(x[ , categorical[i]]))[ ,-1])

                        } # adding the dummy variables into the x matrix or data frame

          x <- x[ , -categorical] # remove the initial categorical variables
   }

  }
  
  }
 


  pmt <- read.csv(cspacedatname, header = TRUE)
  # read in the data

  rw <- which(pmt[ , 1] == 1)
  # find the chosen scheme


  tc <- pmt[rw, -1]
  # the chosen scheme


  dpmt <- pmt[ , -1]
  # all the schemes without the chosen indicator variable


  for(i in 1:dim(dpmt)[1]){

  	        for(j in 1:dim(dpmt)[2]){

                        if (dpmt[i, j]==0) {

           	                            dpmt[i, j] <- -1

                                            }

  	                                }

                          }
    # transform 0 in the constrained space into -1



  if (!all(apply(dpmt, 1, sum)==0)) {

    warning("Warning: test might be anti-conservative if number of intervention clusters does not equal number of control clusters.")
    ## If there is any scheme with unequal cluster assignment, the cptest function has inflated test size.

  }



  if(!is.null(z)){
    fm <- as.formula(paste0("outcome~", paste(names(x), collapse = "+")))
  }else{
    fm <- as.formula("outcome ~ 1")
  }
   
   # the formula of the model

   x$outcome <- outcome
   # merge the outcome and the covariates into one data frame

   if (outcometype == "continuous"){

       ADJMeans <- tapply(lm(formula = fm, data = x)$res, clustername, mean)
       # for continuous outcome, we use linear regression

       Diffs <- as.matrix(dpmt) %*% ADJMeans
       # the permutation statistic

       pvalue <- mean(abs(Diffs) >= abs(Diffs[rw, ]))

   } else if (outcometype == "binary") {

       ADJMeans <- tapply(residuals(glm(formula = fm, family = "binomial", data = x), type = "response"), clustername, mean)
       # for binary outcome, we use logistic regression

       Diffs <- as.matrix(dpmt) %*% ADJMeans
       # the permutation statistic

       pvalue <- mean(abs(Diffs) >= abs(Diffs[rw, ]))

   } else {

   	stop("Error: Please specify the correct outcometype for continuous or binary")
    # for other type outcome, just cannot run the cptest function

   }

 FinalScheme <- data.frame(Cluster_ID = names(ADJMeans), Intervention = as.vector(as.matrix(tc)))
 # the final chosen scheme in the constrained space

 pvalue_message <- paste("Clustered permutation test p-value =", round(pvalue, 4), collapse=' ')
 # the p-value from the permutation test

  return(list(FinalScheme = FinalScheme,
              pvalue = round(pvalue, 4),
  	          pvalue_statement = pvalue_message))

  }





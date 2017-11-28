#' Clustered permutation test for cluster randomized trials
#' @param outcome a vector specifying the individual outcome
#' @param id a vector specifying cluster name
#' @param x a matrix or a data frame specifying the individual-level covariates to analyze
#' @param cspacedatname the csv dataset containing the saved randomization space. This dataset contains the permutation matrix, as well as a variable indicating which row of the permutation matrix was saved as the final schemes.
#' @param outcometype the type of regression model that should be run. Options are "continuous" for linear regression fit and "binary" for logistic regression fit.
#' @param categorical a vector specifying categorical (including binary) variables. This can be names of the columns or numbers of indexes of columns)
#' @keywords clustered permutation test, cluster randomized trails
#' @author Hengshi Yu <hengshi@umich.edu>, John A. Gallis <john.gallis@duke.edu>, Fan Li <frank.li@duke.edu>, Elizabeth L. Turner <liz.turner@duke.edu>
#' @description cptest performs clustered permutation test on the individual oucome for cluster
#' randomized trials (CRTs). The type of the outcome can be specified by the user to be "continuous" or
#' "binary".
#'
#' With specified outcome type being "continuous" or "binary", linear regression or logistic regression is applied on the outcome and the covariates specified for all individuals. Cluster residual means are collected. With the constrained space,
#' the contrast statistic is created from the schemes and the cluster residual means. The permutation test is then conducted on the contrast statistic for the scheme actually utilized.
#' @references
#' Gail, M.H., Mark, S.D., Carroll, R.J., Green, S.B. and Pee, D., 1996. On design considerations and randomization based inference for community intervention trials. Statistics in medicine, 15(11), pp.1069-1092.
#'
#' Li, F., Lokhnygina, Y., Murray, D.M., Heagerty, P.J. and DeLong, E.R., 2016. An evaluation of constrained randomization for the design and analysis of group randomized trials. Statistics in medicine, 35(10), pp.1565-1579.
#' @export
#' @examples
#'\dontrun{
#' Analysis_result <- cptest(outcome = Dickinson_outcome$outcome,
#'                           id = Dickinson_outcome$county,
#'                           x = data.frame(Dickinson_outcome[, c(-1, -7)]),
#'                           cspacedatname="dickinson_constrained.csv",
#'                           outcometype="binary",
#'                           categorical=c("location","incomecat"))
#'}
#'
#' @return \code{FinalScheme} the final scheme in the permutation matrix
#' @return \code{pvalue} the p-value of the intervention effect from the clustered permutation test
#' @return \code{pvalue_statement} the statement of the p-value of the intervention effect from the clustered permutation test






cptest <- function(outcome, id, x, cspacedatname, outcometype, categorical = NULL){


    if (is.null(categorical)) { # if there are no categorical variables specified

   if (sum(apply(x, 2, function(x) length(unique(x[!is.na(x)]))) <= 1) >= 1) {
     ## check the variables with only one unique value, i.e. no variation

     warning(cat("Warning: each of columns", c(1 : dim(x)[2])[apply(x, 2, function(x) length(unique(x[!is.na(x)]))) <= 1], "has only less than or equal to one unique value and no variation"))
   }




   if (!setequal(c(1 : dim(x)[2])[apply(x, 2, function(x) length(unique(x[!is.na(x)]))) <= 4], categorical)) {
     ## check the categorical variables have been correctly specified

     warning(cat("Warning: each of columns", c(1 : dim(x)[2])[apply(x, 2, function(x) length(unique(x))) <= 4], "has less than or equal to 4 unique observations! Please check all the categorical variables have been correctly specified!"))


   }


    } else {  # If there are categorical variables specified

    p <- length(categorical)

     if (is.character(categorical)) {  ## columns names based categorical variables

      if (sum(apply(x, 2, function(x) length(unique(x[!is.na(x)]))) <= 1) >= 1) {
         ## check the variables with only one unique value, i.e. no variation

         warning(cat("Warning: each of", colnames(x)[apply(x, 2, function(x) length(unique(x[!is.na(x)]))) <= 1], "has only less than or equal to one unique value and no variation"))
       }


      if (!setequal(colnames(x)[apply(x, 2, function(x) length(unique(x[!is.na(x)]))) <= 4], categorical)) {
        ## check the categorical variables have been correctly specified


        warning(cat("Warning: each of", setdiff(colnames(x)[apply(x, 2, function(x) length(unique(x))) <= 4], categorical), "has less than or equal to 4 unique observations! Please check all the categorical variables have been correctly specified!"))
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

           warning(cat("Warning: each of columns", c(1:dim(x)[2])[apply(x, 2, function(x) length(unique(x[!is.na(x)]))) <= 1], "has only less than or equal to one unique value and no variation"))

          }


         if (!setequal(c(1:dim(x)[2])[apply(x, 2, function(x) length(unique(x[!is.na(x)]))) <= 4], categorical)) {
                           ## check the categorical variables have been correctly specified

       warning(cat("Warning: each of columns",setdiff(c(1:dim(x)[2])[apply(x, 2, function(x) length(unique(x))) <= 4], categorical), "has less than or equal to 4 unique observations! Please check all the categorical variables have been correctly specified!"))


     }



         for (i in 1:p) {

               x <- data.frame(x, model.matrix(~factor(x[ , categorical[i]]))[ ,-1])

                        } # adding the dummy variables into the x matrix or data frame

          x <- x[ , -categorical] # remove the initial categorical variables
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

    warning("Error: The constrained space contains schemes with unequal cluster assignment. The cptest function is only valid with equal cluster assignment")
    ## If there is any scheme with unequal cluster assignment, the cptest function is not valid

  }


   fm <- as.formula(paste0("outcome~", paste(names(x), collapse = "+")))
   # the formula of the model

   x$outcome <- outcome
   # merge the outcome and the covariates into one data frame

   if (outcometype == "continuous"){

       ADJMeans <- tapply(lm(formula = fm, data = x)$res, id, mean)
       # for continuous outcome, we use linear regression

       Diffs <- as.matrix(dpmt) %*% ADJMeans
       # the permutation statistic

       pvalue <- mean(abs(Diffs) >= abs(Diffs[rw, ]))

   } else if (outcometype == "binary") {

       ADJMeans <- tapply(glm(formula = fm, family = "binomial", data = x)$residuals, id, mean)
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

 pvalue_message <- paste("Adjusted permutation test p-value =", round(pvalue, 4), collapse=' ')
 # the p-value from the permutation test

  return(list(FinalScheme = FinalScheme,
              pvalue = round(pvalue, 4),
  	          pvalue_statement = pvalue_message))

  }





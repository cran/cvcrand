#' Covariate-constrained randomization for cluster randomized trials
#' @param clustername a vector specifying the identification variable of the cluster. If no cluster identification variable is specified, the default is to label the clusters based on the order in which they appear.
#' @param x a data frame specifying the values of cluster-level covariates to balance. With K covariates and n clusters, it will be dimension of \code{n} by \code{K}.
#' @param categorical a vector specifying categorical (including binary) variables. This can be names of the columns or number indexes of columns, but cannot be both. Suppose there are \code{p} categories for a categorical variable, \code{cvcrand} function creates \code{p-1} dummy variables and drops the reference level if the variable is specified as a factor. Otherwise, the first level in the alphanumerical order will be dropped. The results are sensitive to which level is excluded. If the user wants to specify a different level to drop for a \code{p}-level categorical variable, the user can create \code{p-1} dummy variables and these can instead be supplied as covariates to the \code{cvcrand} function. Then, the user needs to specify the dummy variables created to be \code{categorical} when running \code{cvcrand}. In addition, the user could also set the variable as a factor with the specific reference level. If the \code{weights} option is used, the weights for a categorical variable will be replicated on all the dummy variables created.
#' @param weights a vector of user-specified weights for the covariates to calculate the balance score. The weight for a categorical variable will be replicated for the dummy variables created. Note that the \code{weights} option can be used to conduct stratification on variables. For example, a variable with a relatively large weight like \code{1000} and all other variables with a weight of \code{1} will cause the randomization scheme chosen to be stratified by the variable with the large weight, assuming a low \code{cutoff} value is specified.
#' @param ntotal_cluster the total number of clusters to be randomized. It must be a positive integer and equal to the number of rows of the data.
#' @param ntrt_cluster  the number of clusters that the researcher wants to assign to the treatment arm. It must be a positive integer less than the total number of clusters.
#' @param cutoff quantile cutoff of the distribution of balance score below which a randomization scheme is sampled. Its default is \code{0.1}, and it must be between 0 and 1. The \code{cutoff} option is overridden by the \code{numschemes} option.
#' @param numschemes number of randomization schemes to form the constrained space for the final randomization scheme to be selected. If specified, it overrides the option \code{cutoff} and the program will randomly sample the final randomization scheme from the constrained space of randomization schemes with the \code{numschemes} smallest balance scores. It must be a positive integer.
#' @param size number of randomization schemes to simulate if the number of all possible randomization schemes is over \code{size}. Its default is \code{50,000}, and must be a positive integer. It can be overridden by the \code{nosim} option.
#' @param stratify categorical variables on which to stratify the randomization. It overrides the option \code{weights} when specified. This list of categorical variables should be a subset of the \code{categorical} option if specified.
#' @param seed seed for simulation and random sampling. It is needed so that the randomization can be replicated. Its default is \code{12345}.
#' @param balancemetric balance metric to use. Its choices are \code{"l1"} and \code{"l2"}. The default is \code{"l2"}.
#' @param nosim if \code{TRUE}, it overrides the default procedure of simulating when the number of all possible randomization schemes is over the \code{size}, and the program enumerates all randomization schemes. Note: this may consume a lot of memory and cause R to crash
#' @param savedata saves the data set of the constrained randomization space in a csv file if specified by \code{savedata}. The first column of the csv file is an indicator variable of the final randomization scheme in the constrained space. The constrained randomization space will be needed for analysis after the cluster randomized trial is completed if the clustered permutation test is used.
#' @param bhist if \code{TRUE} of the default value, it produces the histogram of all balance scores with a red line on the graph indicating the selected cutoff.
#' @param check_validity boolean argument to check the randomization validity or not
#' @param samearmhi clusters assigned to the same arm as least this often are displayed. The default is \code{0.75}. 
#' @param samearmlo clusters assigned to the same arm at most this often are displayed. The default is \code{0.25}. 
#' @keywords cluster-randomized-trails covariate-constrained-randomization
#' @author Hengshi Yu <hengshi@umich.edu>, Fan Li <fan.f.li@yale.edu>, John A. Gallis <john.gallis@duke.edu>, Elizabeth L. Turner <liz.turner@duke.edu>
#' @description \code{cvrall} performs constrained randomization for cluster randomized
#' trials (CRTs), especially suited for CRTs with a small number of clusters. In constrained randomization,
#' a randomization scheme is randomly sampled from a subset of all possible randomization schemes
#' based on the value of a balancing criterion called a balance score. The \code{cvrall} function has two choices of "l1" and "l2" metrics for balance score.
#'
#' The \code{cvrall} function enumerates all randomization schemes or chooses the unique ones among some simulated randomization schemes as specified by the user.
#' Some cluster-level continuous or categorical covariates are then used to calculate the balance scores for the unique schemes. A subset of the randomization schemes is chosen based on a user-specified cutoff at a certain quantile of the distribution of the balance scores or based on a fixed number of schemes with the smallest balance scores. The \code{cvrall} function
#' treats the subset as the constrained space of randomization schemes and samples one scheme from the constrained space as the final chosen scheme.
#'
#' @references
#' Raab, G.M. and Butcher, I., 2001. Balance in cluster randomized trials. Statistics in medicine, 20(3), pp.351-365.
#'
#' Li, F., Lokhnygina, Y., Murray, D.M., Heagerty, P.J. and DeLong, E.R., 2016. An evaluation of constrained randomization for the design and analysis of group randomized trials. Statistics in medicine, 35(10), pp.1565-1579.
#'
#' Li, F., Turner, E. L., Heagerty, P. J., Murray, D. M., Vollmer, W. M., & DeLong, E. R. (2017). An evaluation of constrained randomization for the design and analysis of group randomized trials with binary outcomes. Statistics in medicine, 36(24), 3791-3806.
#'
#' Gallis, J.A., Li, F., Yu, H. and Turner, E.L., 2018. cvcrand and cptest: Commands for efficient design and analysis of cluster randomized trials using constrained randomization and permutation tests. The Stata Journal, 18(2), pp.357-378.
#'
#' Dickinson, L. M., Beaty, B., Fox, C., Pace, W., Dickinson, W. P., Emsermann, C., & Kempe, A. (2015). Pragmatic cluster randomized trials using covariate constrained randomization: A method for practice-based research networks (PBRNs). The Journal of the American Board of Family Medicine, 28(5), 663-672.
#'
#' Bailey, R.A. and Rowley, C.A., 1987. Valid randomization. Proceedings of the Royal Society of London. A. Mathematical and Physical Sciences, 410(1838), pp.105-124.
#'
#' @export
#' @examples
#'
#'
#' # cvrall examples
#'
#' Design_result <- cvrall(clustername = Dickinson_design$county,
#'                          balancemetric = "l2",
#'                          x = data.frame(Dickinson_design[ , c("location", "inciis",
#'                               "uptodateonimmunizations", "hispanic", "incomecat")]),
#'                          ntotal_cluster = 16,
#'                          ntrt_cluster = 8,
#'                          categorical = c("location", "incomecat"),
#'                          ###### Option to save the constrained space ######
#'                          # savedata = "dickinson_constrained.csv",
#'                          bhist = TRUE,
#'                          cutoff = 0.1,
#'                          seed = 12345, 
#'                          check_validity = TRUE)
#'
#' # cvrall example with weights specified
#'
#' Design_result <- cvrall(clustername = Dickinson_design$county,
#'                          balancemetric = "l2",
#'                          x = data.frame(Dickinson_design[ , c("location", "inciis",
#'                              "uptodateonimmunizations", "hispanic", "incomecat")]),
#'                          ntotal_cluster = 16,
#'                          ntrt_cluster = 8,
#'                          categorical = c("location", "incomecat"),
#'                          weights = c(1, 1, 1, 1, 1),
#'                          cutoff = 0.1,
#'                          seed = 12345, 
#'                          check_validity = TRUE)
#'
#' # Stratification on location, with constrained
#' # randomization on other specified covariates
#'
#'  Design_stratified_result <- cvrall(clustername = Dickinson_design$county,
#'                                      balancemetric = "l2",
#'                                      x = data.frame(Dickinson_design[ , c("location", "inciis",
#'                                          "uptodateonimmunizations", "hispanic", "incomecat")]),
#'                                      ntotal_cluster = 16,
#'                                      ntrt_cluster = 8,
#'                                      categorical = c("location", "incomecat"),
#'                                      weights = c(1000, 1, 1, 1, 1),
#'                                      cutoff = 0.1,
#'                                      seed = 12345)
#'
#'  # An alternative and equivalent way to stratify on location
#'
#'  Design_stratified_result <- cvrall(clustername = Dickinson_design$county,
#'                                      balancemetric = "l2",
#'                                      x = data.frame(Dickinson_design[ , c("location", "inciis",
#'                                          "uptodateonimmunizations", "hispanic", "incomecat")]),
#'                                      ntotal_cluster = 16,
#'                                      ntrt_cluster = 8,
#'                                      categorical = c("location", "incomecat"),
#'                                      stratify = "location",
#'                                      cutoff = 0.1,
#'                                      seed = 12345)
#'
#'  # Stratification on income category
#'  # Two of the income categories contain an odd number of clusters
#'  # Stratification is not strictly possible
#'
#'  Design_stratified_inc_result <- cvrall(clustername = Dickinson_design$county,
#'                                          balancemetric = "l2",
#'                                          x = data.frame(Dickinson_design[ , c("location", "inciis",
#'                                              "uptodateonimmunizations", "hispanic", "incomecat")]),
#'                                          ntotal_cluster = 16,
#'                                          ntrt_cluster = 8,
#'                                          categorical = c("location", "incomecat"),
#'                                          stratify = "incomecat",
#'                                          cutoff = 0.1,
#'                                          seed = 12345)
#' 
#'
#' @return \code{balancemetric} the balance metric used
#' @return \code{allocation} the allocation scheme from constrained randomization
#' @return \code{bscores} the histogram of the balance score with respect to the balance metric
#' @return \code{assignment_message} the statement about how many clusters to be randomized to the intervention and the control arms respectively
#' @return \code{scheme_message} the statement about how to get the whole randomization space to use in constrained randomization
#' @return \code{cutoff_message} the statement about the cutoff in the constrained space
#' @return \code{choice_message} the statement about the selected scheme from constrained randomization
#' @return \code{data_CR} the data frame containing the allocation scheme, the \code{clustername}, and the original data frame of covariates
#' @return \code{baseline_table} the descriptive statistics for all the variables by the two arms from the selected scheme
#' @return \code{cluster_coincidence} cluster coincidence matrix
#' @return \code{cluster_coin_des} cluster coincidence descriptive
#' @return \code{clusters_always_pair} pairs of clusters always allocated to the same arm.
#' @return \code{clusters_always_not_pair} pairs of clusters always allocated to different arms.
#' @return \code{clusters_high_pair}  pairs of clusters randomized to the same arm at least \code{samearmhi} of the time.
#' @return \code{clusters_low_pair} pairs of clusters randomized to the same arm at most \code{samearmlo} of the time.
#' @return \code{overall_allocations} frequency of acceptable overall allocations.


cvrall = function(clustername = NULL, x, categorical = NULL, weights = NULL, ntotal_cluster, ntrt_cluster, cutoff = 0.1, numschemes = NULL, size = 50000, stratify = NULL, seed = NULL, balancemetric = "l2", nosim = FALSE, savedata = NULL, bhist = TRUE, check_validity = FALSE, samearmhi = 0.75, samearmlo = 0.25){

  if (!is.null(seed)) {
      set.seed(seed)
      } else {
      set.seed(12345)
      } ##default seed to be 12345

  if(is.null(clustername)){

    clustername <- 1:dim(x)[1]

  }

  nsub <- length(clustername) # number of clusters

  if (ntotal_cluster %% 1 != 0) {
    stop("Error: number of clusters specified must be an integer.")
  }

  if (ntotal_cluster <= 0) {
    stop("Error: number of clusters specified must be a positive number.")
  }

  if (ntrt_cluster %% 1 != 0) {
    stop("Error: number of treatment clusters specified must be an integer.")
  }

  if (ntrt_cluster <= 0) {
    stop("Error: number of treatment clusters specified must be a positive number.")
  }

   if(nsub != dim(x)[1]) {
    stop("Error: the clustername vector and the covariates data frame must match in the number of clusters")
  }

  if(nsub != ntotal_cluster) {
    stop("Error: number of clusters specified must equal the number of rows in the data set.")
  }

  if (ntotal_cluster <= ntrt_cluster) {
    stop("Error: number of clusters in the treatment arm must be smaller than the total number of clusters!")
  }

  if (cutoff <= 0 | cutoff >= 1) {
    stop("Error: cutoff must be greater than 0 and smaller than 1")
  }

  if ((!is.null(stratify))){

  if (is.null(categorical)){

    stop("Error: categorical must be specified when the stratify option is specified")

  } else if (!all(stratify %in% categorical)){

    stop("Error: the categorical variables for stratification must be a subset of the categorical option")

  }

  }

  if(!is.null(stratify)){

    weights <- rep(1, dim(x)[2]) # give weight of 1 to the variables other than the stratification variables

    if (is.character(stratify)) {

      weights[colnames(x) %in% stratify] <- 1000 # give weight of 1000 to each of the stratification variables

     } else if (is.numeric(stratify)) {

      weights[x] <- 1000

       }

    }

  data <- x ## store the original data

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



      if (!is.null(weights)) { ## if there is user-specified weight vector

              for(i in 1:p){

            x <- data.frame(x, model.matrix(~factor(x[ , categorical[i]]))[ , -1])
                            # add the dummy variables into the x matrix

            weights <- c(weights, rep(weights[which(categorical[i]==colnames(x))], dim(model.matrix(~factor(x[ , categorical[i]])))[2] - 1))
                           # set the weights of the dummy variables to be the same as the correpsonding original categorical variables'
              }

            weights <- weights[-which(colnames(x) %in% categorical)] # removing the initial categorical variables' weights

         } else {   ## if there is no user-specified weights vector

              for(i in 1:p){

             x <- data.frame(x, model.matrix(~factor(x[ , categorical[i]]))[ , -1])
                           # add the dummy variables into the x matrix or data frame
                            }
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

       warning(cat("Warning: each of columns",setdiff(c(1:dim(x)[2])[apply(x, 2, function(x) length(unique(x))) <= 4], categorical), "has less than or equal to 4 unique observations. Please check all the categorical variables have been correctly specified."))


     }

     if (!is.null(weights)) {

       for (i in 1:p) {

         x <- data.frame(x, model.matrix(~factor(x[, categorical[i]]))[,-1])
         # adding the dummy variables into the x matrix or data frame

         weights <- c(weights,rep(weights[categorical[i]], dim(model.matrix(~factor(x[ , categorical[i]])))[2] - 1))
       }
       weights <- weights[-categorical] # removing the initial categorical variable weights

     } else {

       for (i in 1:p) {

         x <- data.frame(x, model.matrix(~factor(x[ , categorical[i]]))[ ,-1])

       } # adding the dummy variables into the x matrix or data frame
     }
          x <- x[ , -categorical] # remove the initial categorical variables
   }

  }

      x<-as.matrix(x)

      np <- dim(x)[2]      # number of covariates for constrained randomization


  id = clustername
  n = ntotal_cluster
  ntrt = ntrt_cluster

   if (choose(nsub, ntrt_cluster) <= size | nosim == TRUE) {       # enumeration if there are not too many clusters

    sim <- 0                     # indicate enumeration
    emr <- t(combn(nsub, ntrt_cluster))      # all the schemes
    R <- dim(emr)[1]
    pmt <- matrix(NA, R, nsub)         # indicating clusters to get the treatment
    S <- R

    for (r in 1:R) {
      pmt[r, emr[r, ]] <- 1
      pmt[r, -emr[r, ]] <- 0
    }

    # calculating the B score in an equivalent method
    if (is.null(weights)) {

      if (balancemetric=="l1") {

          BL <- apply(abs((pmt %*% ((x - matrix(apply(x, 2, mean), nsub, np, byrow = TRUE)) / matrix(apply(x, 2, sd), nsub, np, byrow = TRUE)))), 1, sum)
          # l1 norm B score: if there are no weights specified, the default weights are the standard deviations

      } else if (balancemetric == "l2") {
        BL <- apply(((pmt %*% ((x - matrix(apply(x, 2, mean), nsub, np, byrow = TRUE)) / matrix(apply(x, 2, sd), nsub, np, byrow = TRUE))))^2, 1, sum)
      # l2 norm B score: if there are no weights specified, the default weights are the variances

      }

      } else {

      if (balancemetric == "l1") {

       BL <- apply(abs((pmt %*% ((x - matrix(apply(x, 2, mean), nsub, np, byrow = TRUE)) * matrix(weights, nsub, np, byrow = TRUE) / matrix(apply(x, 2, sd), nsub, np, byrow = TRUE)))), 1, sum)
      # l1 norm B score: if there are weights specified

      } else if (balancemetric=="l2") {
        BL <- apply(((pmt %*% ((x - matrix(apply(x, 2, mean), nsub, np, byrow = TRUE)) * matrix(weights, nsub, np, byrow = TRUE) / matrix(apply(x, 2, sd), nsub, np, byrow = TRUE))))^2, 1, sum)
      # l2 norm B score: if there are weights specified
      }
      }



      if (is.null(numschemes)) {

      cho <- round(R * cutoff)
        # the cutoff thresholding

      subid <- order(BL)[1:cho]
      #  B score subset


      BLcut <- BL[order(BL)[cho]]
         # B score subset bechmark
      if (cho <= 1) {
        ## check the variables with only one unique value, i.e. no variation

        warning("Warning: less than or equal to 1 scheme meets the constraints for the constrained space")
        BLcut = NaN
      }

      rw <- sample(subid, 1)
            # choice from B

      BLchoice <- BL[rw]
           # selected scheme's B score

      inter <- pmt[rw, ]
           # allocation from the choice from  B





     } else {

      if (numschemes >= R) {

        stop("Error: the fixed number of schemes in the constrained space for constrained randomization must be smaller or equal to the total number of randomization schemes.")
      }

      cho <- numschemes
      # the fixed-number of schemes thresholding

      subid <- order(BL)[1:cho]
      # B score subset

      BLcut <- BL[order(BL)[cho]]
      # B score subset bechmark
      if (cho <= 1) {
        ## check the variables with only one unique value, i.e. no variation

        warning("Warning: less than or equal to 1 scheme meets the constraints for the constrained space")
        BLcut = NaN
      }

      rw <- sample(subid, 1)
           # choice from B

      BLchoice <- BL[rw]
         # selected scheme's B score

      inter <- pmt[rw, ]
      # allocation from the choice from  B


    }


    qmt <- pmt[subid, ,drop = FALSE]
    R_result <- dim(qmt)[1]
    summary_constraints <- as.data.frame(cbind(S, R, R_result, paste0(round(R_result/R, 4) * 100, "%")))
    colnames(summary_constraints) <- c("overall allocations", "checked allocations", "accepted allocations", "overall % acceptable")


    ## histogram of Bscores mean sd, B score benchmark, chosen B score, min max  p5 p10 p20 p25 p30 p50 p75 p95

    BL_Quantiles <- quantile(BL, prob = c(0, 0.05, 0.1, 0.2, 0.25, 0.3, 0.5, 0.75, 0.95, 1))

    BL_Quantiles <- c(BLchoice, BLcut, mean(BL), sd(BL), BL_Quantiles)

    BL_Quantiles <- round(BL_Quantiles, 3)

    names(BL_Quantiles)[c(1, 2, 3, 4, 5, 14)] <- c("score (selected scheme)", "cutoff score", "Mean", "SD", "Min", "Max")



    } else {

    sim <- 1
      # indicate simulation

    S <- size
     # randomization sample size

    pmt <- matrix(NA, S, nsub)
    # indicating clusters to get the treatment

    for (s in 1:S) {

      trt <- sample(1:nsub, ntrt_cluster)

      pmt[s, trt] <- 1

      pmt[s, -trt] <- 0
    }

    pmt <- unique(pmt)
     # indicator matrix

    R <- dim(pmt)[1]
    ## number of unique schemes from simulation


    # calculating the B score in an equivalent method
    if (is.null(weights)) {

      if (balancemetric=="l1") {

        BL <- apply(abs((pmt %*% ((x - matrix(apply(x, 2, mean), nsub, np, byrow = TRUE)) / matrix(apply(x, 2, sd), nsub, np, byrow = TRUE)))), 1, sum)
        # l1 norm B score: if there are no weights specified, the default weights are the standard deviations

      } else if (balancemetric=="l2") {BL <- apply(((pmt %*% ((x - matrix(apply(x, 2, mean), nsub, np, byrow = TRUE)) / matrix(apply(x, 2, sd), nsub, np, byrow = TRUE))))^2, 1, sum)
      # l2 norm B score: if there are no weights specified, the default weights are the variances
      }

      } else {

      if (balancemetric=="l1") {

        BL <- apply(abs((pmt %*% ((x - matrix(apply(x, 2, mean), nsub, np, byrow = TRUE)) * matrix(weights, nsub, np, byrow = TRUE) / matrix(apply(x, 2, sd), nsub, np, byrow = TRUE)))), 1, sum)
         # l1 norm B score: if there are weights specified

      } else if (balancemetric=="l2") {BL <- apply(((pmt %*% ((x - matrix(apply(x, 2, mean), nsub, np, byrow = TRUE)) * matrix(weights, nsub, np, byrow = TRUE) / matrix(apply(x, 2, sd), nsub, np, byrow = TRUE))))^2, 1, sum)
      # l2 norm B score: if there are weights specified
      }
      }




       if (is.null(numschemes)) {

              cho <- round(R * cutoff)
                # the cutoff thresholding

              subid <- order(BL)[1:cho]
              #  B score subset

              BLcut <- BL[order(BL)[cho]]
              # B score subset bechmark
              if (cho <= 1) {
                ## check the variables with only one unique value, i.e. no variation

                warning("Warning: less than or equal to 1 scheme meets the constraints for the constrained space")
                BLcut = NaN
              }
              rw <- sample(subid, 1)
                 # choice from B

              BLchoice <- BL[rw]
               # selected scheme's B score

              inter <- pmt[rw, ]
               # allocation from the choice from  B



         } else {


      if (numschemes >= R) {

        stop("Error: the fixed number of schemes in the constrained space for constrained randomization must be smaller or equal to the total number of randomization schemes.")
      }

      cho <- numschemes
      # the fixed number of schemes thresholding

      subid <- order(BL)[1:cho]
      # B score subset

      BLcut <- BL[order(BL)[cho]]
       # B score subset bechmark
      if (cho <= 1) {
        ## check the variables with only one unique value, i.e. no variation

        warning("Warning: less than or equal to 1 scheme meets the constraints for the constrained space")
        BLcut = NaN
      }
      rw <- sample(subid, 1)
       # choice from B

      BLchoice <- BL[rw]
        # selected scheme's B score

      inter <- pmt[rw,]
      # allocation from the choice from  B
    }

    qmt <- pmt[subid, ,drop = FALSE]
    R_result <- dim(qmt)[1]
    summary_constraints <- as.data.frame(cbind(S, R, R_result, paste0(round(R_result/R, 4) * 100, "%")))
    colnames(summary_constraints) <- c("overall allocations", "checked allocations", "accepted allocations", "overall % acceptable")


    ## histogram of Bscores mean sd, B score benchmark, chosen B score, min max  p5 p10 p20 p25 p30 p50 p75 p95
    BL_Quantiles <- quantile(BL, prob = c(0, 0.05, 0.1, 0.2, 0.25, 0.3, 0.5, 0.75, 0.95, 1))

    BL_Quantiles <- c(BLchoice, BLcut, mean(BL), sd(BL), BL_Quantiles)

    BL_Quantiles <- round(BL_Quantiles, 3)


    names(BL_Quantiles)[c(1, 2, 3, 4, 5, 14)] <- c("score (selected scheme)", "cutoff score", "Mean", "SD", "Min", "Max")

    }



  if (!is.null(savedata)) {

    SchemeChosen <- rep(0, dim(pmt)[1])

    SchemeChosen[rw] <- 1

    pmt<-cbind(SchemeChosen, pmt)

    write.csv(pmt[subid, ], file = savedata, row.names=FALSE)
  }
      # output the schemes' matrix of pmt


  if (bhist) {


    if (sim == 1) {

      par(mar = c(4, 3, 1.5, 1))


      hist(BL, main = paste("Histogram of balance scores across", R, "unique schemes", "in", S, "simulated schemes"), xlab = NULL, ylab = NULL, breaks=50)

      title(xlab = "Balance score", line = 2)

      title(ylab = "Frequency", line = 2)

          } else {



      par(mar = c(4, 3, 1.5, 1))

      hist(BL, main = paste("Histogram of balance scores across all", R, "schemes"), xlab = NULL, ylab = NULL, breaks=50)

      title(xlab = "Balance score", line = 2)

      title(ylab = "Frequency", line = 2)

      }



    abline(v = BLcut, col = "red")



   if (is.null(numschemes)) {

    mtext(paste("The", cutoff, "quantile of the", balancemetric, "balance score is", round(BLcut, 3)), col = 'black', side = 1, outer = FALSE, line = 3)

       } else {

     mtext(paste("The", numschemes, "schemes", balancemetric, "balance score is", round(BLcut, 3)), col = "black", side = 1, outer = FALSE, line = 3)


      }



    box()


  }
  # output the balance scores across entire randomization space as well as output its histogram


  coin_matrix <- coin_descri <- al_clusters <- alno_clusters  <- hi_clusters <- lo_clusters <- NULL
  if (check_validity) {
      
      n_pair <- t(combn(nsub, 2))      # all the schemes

      same_arm_count <- same_arm_frac <- diff_arm_count <- diff_arm_frac <- rep(NA, dim(n_pair)[1])
      
      for (j in 1:(nsub - 1)) {
        for (k in (j + 1):nsub) {
          same_arm <- sum((qmt[,j] - qmt[, k]) == 0)
          diff_arm <-  R_result - same_arm
          same_prop <- round(same_arm/R_result, 4)
          diff_prop <- round(1.0 - same_prop, 4)
          index_jk <- which(unlist(lapply(1:dim(n_pair)[1], function(i) setequal(c(j,k), n_pair[i, ]))))
          
          same_arm_count[index_jk] <- same_arm
          same_arm_frac[index_jk] <- same_prop
          diff_arm_count[index_jk] <- diff_arm
          diff_arm_frac[index_jk] <- diff_prop
        }
      }
      first_cluster <- id[n_pair[,1]]
      second_cluster <- id[n_pair[,2]]
      coin_matrix <- as.data.frame(cbind(first_cluster, 
                                  second_cluster, 
                                  same_arm_count, 
                                  paste0(same_arm_frac * 100, "%"), 
                                  diff_arm_count, 
                                  paste0(diff_arm_frac * 100, "%")))
      
      
      # Always togerther
      if (sum(same_arm_frac == 1) > 0) {
        alto_index <- which(same_arm_frac == 1)
        alto_clt_pair <- c()
        for (t in 1:length(alto_index)) {
          clt_index <- alto_index[t]
          alto_clts <- paste0(first_cluster[clt_index], " and ", second_cluster[clt_index])
          alto_clt_pair <- c(alto_clt_pair, alto_clts)
        }

        alto_all_pt <- rep("100.0%", length(alto_index))
        
        al_clusters <- as.data.frame(cbind(alto_clt_pair, alto_all_pt))
        colnames(al_clusters) <- c("cluter pair", "% allocs in the same arm")
      }

      # Always not together
      if (sum(same_arm_frac == 0) > 0) {
        alnoto_index <- which(same_arm_frac == 0)
        alnoto_clt_pair <- c()
        for (t in 1:length(alnoto_index)) {
          clt_index <- alnoto_index[t]
          alnoto_clts <- paste0(first_cluster[clt_index], " and ", second_cluster[clt_index])
          alnoto_clt_pair <- c(alnoto_clt_pair, alnoto_clts)
        }

        alnoto_all_pt <- rep("0.0%", length(alnoto_index))
        
        alno_clusters <- as.data.frame(cbind(alnoto_clt_pair, alnoto_all_pt))
        colnames(alno_clusters) <- c("cluter pair", "% allocs in the same arm")
      }
      
      
      # user specify upper bound
      if (sum(same_arm_frac >= samearmhi) > 0) {
        hi_index <- which(same_arm_frac >= samearmhi)
        hi_pair <- c()
        for (t in 1:length(hi_index)) {
          clt_index <- hi_index[t]
          hi_cluts <- paste0(first_cluster[clt_index], " and ", second_cluster[clt_index])
          hi_pair <- c(hi_pair, hi_cluts)
        }

        hi_pt <- paste0(same_arm_frac[hi_index] * 100, "%")
        hi_num <- same_arm_count[hi_index]
        hi_non_num <- diff_arm_count[hi_index]
        hi_clusters <- as.data.frame(cbind(hi_pair, hi_pt, hi_num, hi_non_num))
        colnames(hi_clusters) <- c("cluster pair", "% allocs in same arm", "# allocs in same arm", "# allocs in different arms")
      }
      
      # user specified lower bound
      if (sum(same_arm_frac <= samearmlo) > 0) {
        lo_index <- which(same_arm_frac <= samearmlo)
        lo_pair <- c()
        for (t in 1:length(lo_index)) {
          clt_index <- lo_index[t]
          lo_cluts <- paste0(first_cluster[clt_index], " and ", second_cluster[clt_index])
          lo_pair <- c(lo_pair, lo_cluts)
        }

        lo_pt <- paste0(same_arm_frac[lo_index] * 100, "%")
        lo_num <- same_arm_count[lo_index]
        lo_non_num <- diff_arm_count[lo_index]
        lo_clusters <- as.data.frame(cbind(lo_pair, lo_pt, lo_num, lo_non_num))
        colnames(lo_clusters) <- c("cluster pair", "% allocs in same arm", "# allocs in same arm", "# allocs in different arms")
      }
      
          
      coin_descri <- rbind(
        c(mean(same_arm_count), sd(same_arm_count), 
        quantile(same_arm_count, prob = c(0, 0.25, 0.5, 0.75, 1))
        ), 
        c(mean(same_arm_frac), sd(same_arm_frac), 
        quantile(same_arm_frac, prob = c(0, 0.25, 0.5, 0.75, 1))
        ), 
        c(mean(diff_arm_count), sd(diff_arm_count), 
        quantile(diff_arm_count, prob = c(0, 0.25, 0.5, 0.75, 1))
        ), 
        c(mean(diff_arm_frac), sd(diff_arm_frac), 
        quantile(diff_arm_frac, prob = c(0, 0.25, 0.5, 0.75, 1))
        )
      )
      coin_descri <- round(coin_descri, 3)

      


      colnames(coin_descri) <- c("Mean", "Std Dev", "Minimum", "25th Pctl", 
                                 "Median", "75th Pctl", "Maximum")

      rownames(coin_descri) <- c("samecount", "samefrac", "diffcount", "difffrac")


      colnames(coin_matrix) <- c("cluster 1", 
                                        "cluster 2", 
                                        "# same arm", 
                                        "% same arm", 
                                        "# different arms", 
                                        "% different arms")
      
    }

    # the allocation of schemes from l1 norm B score or l2 norm B score
    allocation <- cbind(clustername, inter)

    colnames(allocation) <- c("clustername", "allocation")

    allocation <- as.data.frame(allocation)

    BL_Quantiles <- data.frame(names(BL_Quantiles), BL_Quantiles)

    rownames(BL_Quantiles) <- NULL

    colnames(BL_Quantiles) <- NULL

    assignment_message <- paste("You have indicated that you want to assign", ntrt_cluster, "clusters to treatment", "and", ntotal_cluster - ntrt_cluster, "to control")

    # indicate a enumeration process or a simulation process with the detailed number of schemes
    if (sim == 1) {

             scheme_message<-paste("Simulating", S, "schemes with", R, "unique schemes for", ntrt_cluster, "clusters in the treatment arm out of", ntotal_cluster, "clusters in total")

          } else {

    scheme_message<-paste("Enumerating all the", R, "schemes for", ntrt_cluster, "clusters in the treatment arm out of", ntotal_cluster, "clusters in total")

      }


      # indicate the cutoff for constrained randomization
     if (is.null(numschemes)) {

       cutoff_message <- paste("The quantile cutoff value is", cutoff, "based on the", balancemetric, "balance metric, the cutoff balance score is", round(BLcut, 3))

       } else {

       cutoff_message <- paste("The fixed number of schemes is", numschemes, "based on the", balancemetric, "balance metric, the cutoff balance score is", round(BLcut, 3))
      }


      # indicate the chosen scheme's BL from the BL criterion
       choice_message <- paste("Balance score of selected scheme by", balancemetric ,"is", round(BLchoice, 3), collapse=' ')



     data_merge <- data.frame(inter, clustername, data)
      # data frame including the chosen scheme from BL, cluster name and the covariates

     colnames(data_merge)[1] <- "arm"

     if (!is.null(categorical)) { 
     # put the categorical variables into factors to prepare for the "CreateTableOne" function

        if (is.character(categorical)) {
                   varsToFactor <- categorical
         # if the categorical variables' names have been specified, these are the names in the data frame for the variables to be processed to be factors.

       } else {

        varsToFactor <- colnames(data_merge)[categorical + 2]
      # if the categorical variables' indexes in the original covariates' matrix or data frame have been specified, the names in the data frame can be extracted to be processed to be factors.
     }


      data_merge[varsToFactor] <- lapply(data_merge[varsToFactor], factor)
      # put the categorical variables to factors

    }

      Descriptive_statistics <- CreateTableOne( vars = colnames(data_merge)[c(-1, -2)], strata=c("arm"), data=data_merge, test =  FALSE)
       # create the descriptive table to compare the two arms from constrained randomization

     invisible(capture.output(DS <- as.data.frame(print(Descriptive_statistics))))
     # transform the table into a data frame

     colnames(DS)[1:2] <- c("arm = 0", "arm = 1")

  return(list(balancemetric = balancemetric,
              allocation = allocation,
              bscores = BL_Quantiles,
              assignment_message = assignment_message,
              scheme_message = scheme_message,
              cutoff_message = cutoff_message,
              choice_message = choice_message,
              data_CR = data_merge,
              baseline_table = DS, 
              cluster_coincidence = coin_matrix, 
              cluster_coin_des = coin_descri, 
              clusters_always_pair = al_clusters, 
              clusters_always_not_pair = alno_clusters, 
              clusters_high_pair = hi_clusters, 
              clusters_low_pair = lo_clusters, 
              overall_allocations = summary_constraints))

  ## return the allocations, number of schemes, cutoff and choice messages, resulted arms comparisons from BL

  }



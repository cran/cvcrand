#' Constrained randomization for cluster randomized trials
#' @param id a vector specifying cluster name
#' @param x a data frame specifying the cluster-level covariates to balance
#' @param categorical a vector specifying categorical (including binary) variables. This can be names of the columns or number indexes of columns. Suppose there are \code{K} categories for a categorical variable, \code{cvcrand} function creates \code{K-1} dummy variables and drops the first level based on the alphanumeric order. If the user wants to specify different level to drop for a categorical variable, the user can create \code{K-1} dummy variables as wished to use before running \code{cvcrand}
#' @param weight a vector of user-specified weights for the covariates
#' @param n the total number of clusters. It must be a positive integer
#' @param ntrt the number of clusters in the treatment arm; must be a positive integer less than \code{n}
#' @param cutoff cutoff of distribution of balance scores below which a randomization scheme is sampled. Its default is 0.1, and must be between 0 and 1.
#' @param numschemes number of schemes to sample. If specified, it overrides cutoff and it must be a positive integer.
#' @param size number of simulated schemes when simulation is invoked. Its default is 50,000, and must be a positive integer.
#' @param seed randomization seed. Its default is 12345
#' @param metric balance metric to use. Its choices are "L1" and "L2". The default is "L2"
#' @param nosim If TRUE, it makes the program enumerate all randomization schemes, regardless of number
#' @param savedata the path of the csv file of the constrained randomization space if specified. The first column of the csv data is an indicator variable of the selected randomization scheme. The csv file is needed for clustered permutation test analysis
#' @param savebscores the path of the csv file of the vector of balance scores across the entire randomization space if specified. When this option is specified, a histogram is also produced which displays the distribution of all balance scores with a red line on the graph indicating the selected cutoff
#' @keywords cluster randomized trails, constrained randomization
#' @author Hengshi Yu <hengshi@umich.edu>, John A. Gallis <john.gallis@duke.edu>, Fan Li <frank.li@duke.edu>, Elizabeth L. Turner <liz.turner@duke.edu>
#' @description cvcrand performs constrained randomization for cluster randomized
#' trials (CRTs), especially with a small number of clusters. In constrained randomization,
#' a randomized scheme is randomly sampled from a subset of all possible randomization schemes
#' based on the value of a balancing criterion proposed by Raab and Butcher (2001) called balance score. And \code{cvcrand} function has two choices of "L1" and "L2" metrics for balance score.
#'
#' \code{cvcrand} function enumerates all randomization schemes or chooses the unique ones among some simulated randomization schemes as specified by the user.
#' Some cluster-level "continuous" or "categorical" covariates are then used to calculate the balance scores for the focused schemes. A subset of the randomization schemes is chosen based on user-specified cutoff of certain quantile of balance scores or some fixed number of schemes with the smallest balance scores. \code{cvcrand} function
#' treats the subset as the constrained space of randomization schemes and samples one scheme from the constrained space as the final chosen scheme.
#'
#' @references
#' Raab, G.M. and Butcher, I., 2001. Balance in cluster randomized trials. Statistics in medicine, 20(3), pp.351-365.
#'
#' Li, F., Lokhnygina, Y., Murray, D.M., Heagerty, P.J. and DeLong, E.R., 2016. An evaluation of constrained randomization for the design and analysis of group randomized trials. Statistics in medicine, 35(10), pp.1565-1579.
#' @export
#' @examples
#' Design_result<-cvcrand(id = Dickinson_design$county,
#'                metric = "L2",
#'                x = data.frame(Dickinson_design[, -1]),
#'                n = 16,
#'                ntrt = 8,
#'                categorical = c("location", "incomecat"),
#'                savedata = "dickinson_constrained.csv",
#'                savebscores = "dickinson_bscores.csv",
#'                cutoff = 0.1,
#'                seed = 12345)
#'
#' @return \code{metric} the metric we use
#' @return \code{Allocation} the allocation scheme from the constrained randomization
#' @return \code{Bscores} the histogram of the B score with respect to the metric
#' @return \code{assignment_message} the statement about how many clusters to randomize to the intervention and the control arms respectively
#' @return \code{scheme_message} the statement about how to get the whole randomization space of schemes
#' @return \code{BL_cut_message} the statement about the benchmark of the constrained space
#' @return \code{BLchoice_message} the statement about the selected scheme from constrained randomization
#' @return \code{data_BL} the data frame containing the allocation scheme, the id as well as the original covariates' matrix or data frame
#' @return \code{BL_result} the descriptive statistics for all the variables in the original covariates' matrix or data frame in the two arms from the selected scheme

cvcrand = function(id = NULL, x, categorical = NULL, weight = NULL, n, ntrt, cutoff = 0.1, numschemes = NULL, size = 50000, seed = NULL, metric = "L2", nosim = FALSE, savedata = NULL, savebscores = NULL){

  if (!is.null(seed)) {
      set.seed(seed)
      } else {
      set.seed(12345)
      } ##default seed to be 12345

  if(is.null(id)){

    id <- 1:dim(x)[1]

  }

  nsub <- length(id) # number of clusters

  if (n %% 1 != 0) {
    stop("Error: number of clusters is not an integer!")
  }

  if (n <= 0) {
    stop("Error: number of clusters is not a positive number!")
  }

  if (ntrt %% 1 != 0) {
    stop("Error: number of treatment clusters is not an integer!")
  }

  if (ntrt <= 0) {
    stop("Error: number of treatment clusters is not a positive number!")
  }

   if(nsub != dim(x)[1]) {
    stop("Error: the id vector and the covariates matrix or data frame do not match in the number of clusters")
  }

  if(nsub != n) {
    stop("Error: number of clusters is not correctly set!")
  }

  if (n <= ntrt) {
    stop("Error: number of clusters in the treatment arm is greater than or equal to the total number of clusters!")
  }

  if (cutoff <= 0 | cutoff >= 1) {
    stop("Error: cutoff should be greater than 0 and smaller than 1")
  }


  data <- x ## store the original data

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


      if (!is.null(weight)) { ## if there is user-specified weight vector

              for(i in 1:p){

            x <- data.frame(x, model.matrix(~factor(x[ , categorical[i]]))[ , -1])
                            # add the dummy variables into the x matrix

            weight <- c(weight, rep(weight[which(categorical[i]==colnames(x))], dim(model.matrix(~factor(x[ , categorical[i]])))[2] - 1))
                           # set the weight of the dummy variables to be the same as the correpsonding original categorical variables'
              }

            weight<-weight[-which(colnames(x) %in% categorical)] # removing the initial categorical variables' weights

         } else {   ## if there is no user-specified weight vector

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

           warning(cat("Warning: each of columns", c(1:dim(x)[2])[apply(x, 2, function(x) length(unique(x[!is.na(x)]))) <= 1], "has only less than or equal to one unique value and no variation"))

          }


         if (!setequal(c(1:dim(x)[2])[apply(x, 2, function(x) length(unique(x[!is.na(x)]))) <= 4], categorical)) {
                           ## check the categorical variables have been correctly specified

       warning(cat("Warning: each of columns",setdiff(c(1:dim(x)[2])[apply(x, 2, function(x) length(unique(x))) <= 4], categorical), "has less than or equal to 4 unique observations! Please check all the categorical variables have been correctly specified!"))


     }

     if (!is.null(weight)) {

       for (i in 1:p) {

         x <- data.frame(x, model.matrix(~factor(x[, categorical[i]]))[,-1])
         # adding the dummy variables into the x matrix or data frame

         weight <- c(weight,rep(weight[categorical[i]], dim(model.matrix(~factor(x[ , categorical[i]])))[2] - 1))
       }
       weight <- weight[-categorical] # removing the initial categorical variable weights

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




   if (choose(nsub, ntrt) <= size | nosim == TRUE) {       # enumeration if there are not too many clusters

    sim <- 0                     # indicate enumeration
    emr <- t(combn(nsub, ntrt))      # all the schemes
    R <- dim(emr)[1]
    pmt <- matrix(NA, R, nsub)         # indicating clusters to get the treatment
    S <- R

    for (r in 1:R){
      pmt[r, emr[r, ]] <- 1
      pmt[r, -emr[r, ]] <- 0
    }

    # calculating the B score in an equivalent method
    if (is.null(weight)) {

      if (metric=="L1") {

          BL <- apply(abs((pmt %*% ((x - matrix(apply(x, 2, mean), nsub, np, byrow = TRUE)) / matrix(apply(x, 2, sd), nsub, np, byrow = TRUE)))), 1, sum)
          # L1 norm B score: if there are no weights specified, the default weights are the standard deviations

      } else if (metric=="L2") {
        BL <- apply(((pmt %*% ((x - matrix(apply(x, 2, mean), nsub, np, byrow = TRUE)) / matrix(apply(x, 2, sd), nsub, np, byrow = TRUE))))^2, 1, sum)
      # L2 norm B score: if there are no weights specified, the default weights are the variances

      }

      } else {

      if (metric=="L1") {

       BL <- apply(abs((pmt %*% ((x - matrix(apply(x, 2, mean), nsub, np, byrow = TRUE)) * matrix(weight, nsub, np, byrow = TRUE) / matrix(apply(x, 2, sd), nsub, np, byrow = TRUE)))), 1, sum)
      # L1 norm B score: if there are weights specified

      } else if (metric=="L2") {
        BL <- apply(((pmt %*% ((x - matrix(apply(x, 2, mean), nsub, np, byrow = TRUE)) * matrix(weight, nsub, np, byrow = TRUE) / matrix(apply(x, 2, sd), nsub, np, byrow = TRUE))))^2, 1, sum)
      # L2 norm B score: if there are weights specified
      }
      }



      if (is.null(numschemes)) {

      cho <- round(R * cutoff)
        # the cutoff thresholding

      subid <- order(BL)[1:cho]
      #  B score subset


      BLcut <- BL[order(BL)[cho]]
         # B score subset bechmark

      rw <- sample(subid, 1)
            # choice from B

      BLchoice <- BL[rw]
           # selected scheme's B score

      inter <- pmt[rw, ]
           # allocation from the choice from  B





     } else {

      if (numschemes >= R) {

        stop("Error: the fixed number of schemes in the constrained space for constrained randomization is larger than the total number of randomization schemes.")
      }

      cho <- numschemes
      # the fixed-number of schemes thresholding

      subid <- order(BL)[1:cho]
      # B score subset

      BLcut <- BL[order(BL)[cho]]
      # B score subset bechmark

      rw <- sample(subid, 1)
           # choice from B

      BLchoice <- BL[rw]
         # selected scheme's B score

      inter <- pmt[rw, ]
      # allocation from the choice from  B


    }


    ## histogram of Bscores mean sd, B score benchmark, chosen B score, min max  p5 p10 p20 p25 p30 p50 p75 p95

    BL_Quantiles <- quantile(BL, prob = c(0, 0.05, 0.1, 0.2, 0.25, 0.3, 0.5, 0.75, 0.95, 1))

    BL_Quantiles <- c(BLcut, BLchoice, mean(BL), sd(BL), BL_Quantiles)

    names(BL_Quantiles)[c(1, 2, 3, 4, 5, 14)] <- c("CR-Benchmark", "Selected Point", "Mean", "SD", "Min", "Max")



    } else {

    sim <- 1
      # indicate simulation

    S <- size
     # randomization sample size

    pmt <- matrix(NA, S, nsub)
    # indicating clusters to get the treatment

    for (s in 1:S) {

      trt <- sample(1:nsub, ntrt)

      pmt[s, trt] <- 1

      pmt[s, -trt] <- 0
    }

    pmt <- unique(pmt)
     # indicator matrix

    R <- dim(pmt)[1]
    ## number of unique schemes from simulation


    # calculating the B score in an equivalent method
    if (is.null(weight)) {

      if (metric=="L1") {

        BL <- apply(abs((pmt %*% ((x - matrix(apply(x, 2, mean), nsub, np, byrow = TRUE)) / matrix(apply(x, 2, sd), nsub, np, byrow = TRUE)))), 1, sum)
        # L1 norm B score: if there are no weights specified, the default weights are the standard deviations

      } else if (metric=="L2") {BL <- apply(((pmt %*% ((x - matrix(apply(x, 2, mean), nsub, np, byrow = TRUE)) / matrix(apply(x, 2, sd), nsub, np, byrow = TRUE))))^2, 1, sum)
      # L2 norm B score: if there are no weights specified, the default weights are the variances
      }

      } else {

      if (metric=="L1") {

        BL <- apply(abs((pmt %*% ((x - matrix(apply(x, 2, mean), nsub, np, byrow = TRUE)) * matrix(weight, nsub, np, byrow = TRUE) / matrix(apply(x, 2, sd), nsub, np, byrow = TRUE)))), 1, sum)
         # L1 norm B score: if there are weights specified

      } else if (metric=="L2") {BL <- apply(((pmt %*% ((x - matrix(apply(x, 2, mean), nsub, np, byrow = TRUE)) * matrix(weight, nsub, np, byrow = TRUE) / matrix(apply(x, 2, sd), nsub, np, byrow = TRUE))))^2, 1, sum)
      # L2 norm B score: if there are weights specified
      }
      }




       if (is.null(numschemes)) {

              cho <- round(R * cutoff)
                # the cutoff thresholding

              subid <- order(BL)[1:cho]
              #  B score subset

              BLcut <- BL[order(BL)[cho]]
              # B score subset bechmark

              rw <- sample(subid, 1)
                 # choice from B

              BLchoice <- BL[rw]
               # selected scheme's B score

              inter <- pmt[rw, ]
               # allocation from the choice from  B



         } else {


      if (numschemes >= R) {

        stop("Error: the fixed number of schemes in the constrained space for constrained randomization is larger than the total number of randomization schemes.")
      }

      cho <- numschemes
      # the fixed number of schemes thresholding

      subid <- order(BL)[1:cho]
      # B score subset

      BLcut <- BL[order(BL)[cho]]
       # B score subset bechmark

      rw <- sample(subid, 1)
       # choice from B

      BLchoice <- BL[rw]
        # selected scheme's B score

      inter <- pmt[rw,]
      # allocation from the choice from  B
    }



    ## histogram of Bscores mean sd, B score benchmark, chosen B score, min max  p5 p10 p20 p25 p30 p50 p75 p95
    BL_Quantiles <- quantile(BL, prob = c(0, 0.05, 0.1, 0.2, 0.25, 0.3, 0.5, 0.75, 0.95, 1))

    BL_Quantiles <- c(BLcut, BLchoice, mean(BL), sd(BL), BL_Quantiles)

    names(BL_Quantiles)[c(1, 2, 3, 4, 5, 14)] <- c("CR-Benchmark", "Selected Point", "Mean", "SD", "Min", "Max")

    }



  if (!is.null(savedata)){

    SchemeChosen <- rep(0, dim(pmt)[1])

    SchemeChosen[rw] <- 1

    pmt<-cbind(SchemeChosen, pmt)

    write.csv(pmt[subid, ], file = savedata, row.names=FALSE)
  }
      # output the schemes' matrix of pmt


  if (!is.null(savebscores)){




    write.csv(BL, file = savebscores, row.names=FALSE)

    if(sim == 1){

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

    mtext(paste(cutoff, "quantile", metric, "balance score benchmark of", round(BLcut, 3)), col = 'red', side = 1, outer = FALSE, line = 3)

       } else {

     mtext(paste(numschemes, "schemes", metric, "balance score benchmark of", round(BLcut, 3)), col = "red", side = 1, outer = FALSE, line = 3)


      }



    box()


  }
  # output the balance scores across entire randomization space as well as output its histogram

    # the allocation of schemes from L1 norm B score or L2 norm B score
    Allocation <- cbind(id, inter)

    colnames(Allocation) <- c("id", "allocation")

    assignment_message <- paste("You have indicated that you want to assign", ntrt, "clusters to treatment", "and", n - ntrt, "to control")

    # indicate a enumeration process or a simulation process with the detailed number of schemes
    if(sim == 1){

             scheme_message<-paste("Simulating", S, "schemes with", R, "unique schemes for", ntrt, "clusters in the treatment arm out of", n, "clusters in total")

          } else {

    scheme_message<-paste("Enumerating all the", R, "schemes for", ntrt, "clusters in the treatment arm out of", n, "clusters in total")

      }


      # indicate the benchmark of BL for constrained randomization
     if (is.null(numschemes)) {

       BLcut_message <- paste("By cutoff quantile of", cutoff, "of", metric, ", the BL constrained randomization benchmark is", round(BLcut, 3))

       } else {

       BLcut_message <- paste("By fixed number of schemes of", numschemes, "of", metric, ", the BL constrained randomization benchmark is", round(BLcut, 3))
      }


      # indicate the chosen scheme's BL from the BL criterion
       BLchoice_message <- paste("Balance score of selected scheme by", metric ,"has B score of", round(BLchoice, 3), collapse=' ')



     data_merge <- data.frame(inter, id, data)
      # data frame including the chosen scheme from BL, cluster id and the covariates


     if(!is.null(categorical)){
     # put the categorical variables into factors to prepare for the "CreateTableOne" function

        if(is.character(categorical)){
                   varsToFactor <- categorical
         # if the categorical variables' names have been specified, these are the names in the data frame for the variables to be processed to be factors.

       } else {

        varsToFactor <- colnames(data_merge)[categorical + 2]
      # if the categorical variables' indexes in the original covariates' matrix or data frame have been specified, the names in the data frame can be extracted to be processed to be factors.
     }


      data_merge[varsToFactor] <- lapply(data_merge[varsToFactor], factor)
      # put the categorical variables to factors

    }

      Descriptive_statistics <- CreateTableOne( vars = colnames(data_merge)[c(-1, -2)], strata=c("inter"), data=data_merge)
       # create the descriptive table to compare the two arms from constrained randomization using BL


  return(list(metric = metric,
             Allocation = Allocation,
             Bscores = BL_Quantiles,
             assignment_message = assignment_message,
             scheme_message = scheme_message,
             BLcut_message = BLcut_message,
             BLchoice_message = BLchoice_message,
             data_BL = data_merge,
             BL_result = Descriptive_statistics))

  ## return the allocations, number of schemes, BL cutoff and choice messages, resulted arms comparisons from BL

  }



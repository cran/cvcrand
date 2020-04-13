#' Covariate-by-covariate constrained randomization for cluster randomized trials
#' @param clustername a vector specifying the identification variable of the cluster. If no cluster identification variable is specified, the default is to label the clusters based on the order in which they appear.
#' @param x a data frame specifying the values of cluster-level covariates to balance. With K covariates and n clusters, it will be dimension of \code{n} by \code{K}.
#' @param categorical a vector specifying categorical (including binary) variables. This can be names of the columns or number indexes of columns, but cannot be both. Suppose there are \code{p} categories for a categorical variable, \code{cvcrand} function creates \code{p-1} dummy variables and drops the reference level if the variable is specified as a factor. Otherwise, the first level in the alphanumerical order will be dropped. The results are sensitive to which level is excluded. If the user wants to specify a different level to drop for a \code{p}-level categorical variable, the user can create \code{p-1} dummy variables and these can instead be supplied as covariates to the \code{cvcrand} function. Then, the user needs to specify the dummy variables created to be \code{categorical} when running \code{cvcrand}. In addition, the user could also set the variable as a factor with the specific reference level. If the \code{weights} option is used, the weights for a categorical variable will be replicated on all the dummy variables created.
#' @param constraints a vector of user-specified constraints for all covariates. \code{"any"} means no constraints. If not \code{"any"}, the first character letter of \code{"m"} denotes absolute mean difference, and \code{"s"} means absolute sum difference. If the second character is \code{"f"}, the previous metric is constrained to be smaller or equal to the fraction with the number followed of the overall mean for \code{"m"} or mean arm total for \code{"s"}. If not \code{"f"} at the second character, the metric is just constrained to be smaller or equal to the value following letter(s). 
#' @param ntotal_cluster the total number of clusters to be randomized. It must be a positive integer and equal to the number of rows of the data.
#' @param ntrt_cluster  the number of clusters that the researcher wants to assign to the treatment arm. It must be a positive integer less than the total number of clusters.
#' @param size number of randomization schemes to simulate if the number of all possible randomization schemes is over \code{size}. Its default is \code{50,000}, and must be a positive integer. It can be overriden by the \code{nosim} option.
#' @param seed seed for simulation and random sampling. It is needed so that the randomization can be replicated. Its default is \code{12345}.
#' @param nosim if TRUE, it overrides the default procedure of simulating when the number of all possible randomization schemes is over \code{size}, and the program enumerates all randomization schemes. Note: this may consume a lot of memory and cause R to crash
#' @param savedata saves the data set of the constrained randomization space in a csv file if specified by \code{savedata}. The first column of the csv file is an indicator variable of the final randomization scheme in the constrained space. The constrained randomization space will be needed for analysis after the cluster randomized trial is completed if the clustered permutation test is used.
#' @param check_validity boolean argument to check the randomization validity or not
#' @param samearmhi clusters assigned to the same arm as least this often are displayed. The default is \code{0.75}. 
#' @param samearmlo clusters assigned to the same arm at most this often are displayed. The default is \code{0.25}. 
#' @keywords cluster-randomized-trails covariate-by-covariate-constrained-randomization
#' @author Hengshi Yu <hengshi@umich.edu>, Fan Li <fan.f.li@yale.edu>, John A. Gallis <john.gallis@duke.edu>, Elizabeth L. Turner <liz.turner@duke.edu>
#' @description \code{cvrcov} performs covariate-by-covariate constrained randomization for cluster randomized
#' trials (CRTs), especially suited for CRTs with a small number of clusters. In constrained randomization,
#' a randomization scheme is randomly sampled from a subset of all possible randomization schemes
#' based on the constraints on each covariate. 
#'
#' The \code{cvrcov} function enumerates all randomization schemes or simulates a fixed size of unique randomization schemes as specified by the user.
#' A subset of the randomization schemes is chosen based on user-specified covariate-by-covariate constraints.  \code{cvrcov} treats the subset as the constrained space 
#' of randomization schemes and samples one scheme from the constrained space as the final chosen scheme.
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
#' Greene, E.J., 2017. A SAS macro for covariate-constrained randomization of general cluster-randomized and unstratified designs. Journal of statistical software, 77(CS1).
#' @export
#' @examples
#'
#' 
#' # cvrcov example
#'
#' Dickinson_design_numeric <- Dickinson_design
#' Dickinson_design_numeric$location = (Dickinson_design$location == "Rural") * 1
#' Design_cov_result <- cvrcov(clustername = Dickinson_design_numeric$county,
#'                             x = data.frame(Dickinson_design_numeric[ , c("location", "inciis",
#'                                 "uptodateonimmunizations", "hispanic", "income")]),
#'                             ntotal_cluster = 16,
#'                             ntrt_cluster = 8,
#'                             constraints = c("s5", "mf.5", "any", "mf0.2", "mf0.2"), 
#'                             categorical = c("location"),
#'                             ###### Option to save the constrained space ######
#'                             # savedata = "dickinson_cov_constrained.csv",
#'                             seed = 12345, 
#'                             check_validity = TRUE)
#' 
#' 
#' @return \code{allocation} the allocation scheme from constrained randomization
#' @return \code{assignment_message} the statement about how many clusters to be randomized to the intervention and the control arms respectively
#' @return \code{scheme_message} the statement about how to get the whole randomization space to use in constrained randomization
#' @return \code{data_CR} the data frame containing the allocation scheme, the \code{clustername}, and the original data frame of covariates
#' @return \code{baseline_table} the descriptive statistics for all the variables by the two arms from the selected scheme
#' @return \code{cluster_coincidence} cluster coincidence matrix
#' @return \code{cluster_coin_des} cluster coincidence descriptive
#' @return \code{clusters_always_pair} pairs of clusters always allocated to the same arm.
#' @return \code{clusters_always_not_pair} pairs of clusters always allocated to different arms.
#' @return \code{clusters_high_pair}  pairs of clusters randomized to the same arm at least \code{samearmhi} of the time.
#' @return \code{clusters_low_pair} pairs of clusters randomized to the same arm at most \code{samearmlo} of the time.
#' @return \code{overall_allocations} frequency of acceptable overall allocations.
#' @return \code{overall_summary} summary of covariates with constraints in the constrained space

cvrcov = function(clustername = NULL, x, categorical = NULL,  constraints, ntotal_cluster, ntrt_cluster,  size = 50000, seed = NULL, nosim = FALSE, savedata = NULL, check_validity = FALSE, samearmhi = 0.75, samearmlo = 0.25){

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



          } else if (is.numeric(categorical)) { ## columns indexes based categorical variables

            if(sum(apply(x, 2, function(x) length(unique(x[!is.na(x)]))) <= 1) >= 1){
             ## check the variables with only one unique value, i.e. no variation

           warning(cat("Warning: each of columns", c(1:dim(x)[2])[apply(x, 2, function(x) length(unique(x[!is.na(x)]))) <= 1], "has only less than or equal to one unique value and no variation."))

          }


         if (!setequal(c(1:dim(x)[2])[apply(x, 2, function(x) length(unique(x[!is.na(x)]))) <= 4], categorical)) {
                           ## check the categorical variables have been correctly specified

       warning(cat("Warning: each of columns",setdiff(c(1:dim(x)[2])[apply(x, 2, function(x) length(unique(x))) <= 4], categorical), "has less than or equal to 4 unique observations. Please check all the categorical variables have been correctly specified."))


     }

  

   }

  }

   x <- as.matrix(x)

      np <- dim(x)[2]      # number of covariates for constrained randomization

  if(length(constraints) != np){

    stop("Error: constraints should be matched with all variables")
  
  }
id = clustername
n = ntotal_cluster
ntrt = ntrt_cluster


if (choose(nsub, ntrt) <= size | nosim == TRUE) {       # enumeration if there are not too many clusters

    sim <- 0                     # indicate enumeration
    emr <- t(combn(nsub, ntrt))      # all the schemes
    R <- dim(emr)[1]
    pmt <- matrix(NA, R, nsub)         # indicating clusters to get the treatment
    qmt <- matrix(NA, R, nsub)
    omt <- matrix(NA, R, nsub)
    mumt <- matrix(NA, R, nsub)
    omumt <- matrix(NA, R, nsub)

    S <- R

    for (r in 1:R){
      omt[r, ] <- 1/2 # for mean arm total

      omumt[r, ] <- 1/nsub  # for overall mean
      
      pmt[r, emr[r, ]] <- 1 # for arm total
      pmt[r, -emr[r, ]] <- -1 # for arm total
      
      mumt[r, emr[r, ]]  <- 1/ntrt # for arm mean
      mumt[r, -emr[r, ]] <- -1/(nsub - ntrt) # for arm mean
      
      
      qmt[r, emr[r, ]] <- 1 # for final
      qmt[r, -emr[r, ]] <- 0 # for final
    }

  any_or_not <- 1 * (constraints == "any")
  constraints_list <- substring(constraints, 1, 1) # "m" or "s"
  fraction_or_not <- substring(constraints, 2, 2)
  constraints_value <- rep(NA, length(constraints))
  constraints_value[fraction_or_not == "f"] = as.numeric(substring(constraints, 3)[fraction_or_not == "f"])  
  constraints_value[fraction_or_not != "f" & fraction_or_not != "n"] = as.numeric(substring(constraints, 2)[fraction_or_not != "f" & fraction_or_not != "n"])

  use_fraction_or_not = 1 * (fraction_or_not == "f")

  sat_num = R
  for(i in 1:length(constraints)){  
    if(any_or_not[i] == 0){
      if(constraints_list[i] == "m"){
        BL_means <- abs(mumt %*% as.matrix(x[,i]))
        if(use_fraction_or_not[i] == 1){ 
          BL_omean <- omumt %*% as.matrix(x[,i])
          sat_list <- c(1 * (BL_means <= constraints_value[i] * BL_omean))
        }else{
          sat_list <- c(1 * (BL_means <= constraints_value[i]))
        }
      }else if(constraints_list[i] == "s"){
        BL_sums <- abs((pmt %*% as.matrix(x[,i])))
        if(use_fraction_or_not[i] == 1){
          BL_osum <- omt %*% as.matrix(x[,i])
          sat_list <- c(1 * (BL_sums <= constraints_value[i] * BL_osum))
        }else{
          sat_list <- c(1 * (BL_sums <= constraints_value[i]))
        }
      }
    }
    sat_num = sum(sat_list)
    if(sat_num > 0){
      omt <- omt[which(sat_list == 1), ]
      omumt <- omumt[which(sat_list == 1), ]
      pmt <- pmt[which(sat_list == 1), ]
      mumt <- mumt[which(sat_list == 1), ]
      qmt <- qmt[which(sat_list == 1), ]
    }else{
      stop("Error: there is not any scheme that satifies the constraints.")
    }   
  }

  if(!is.null(categorical)){
    if(is.character(categorical)){
      cat_index <- which(colnames(x) %in% categorical)
    }else if(is.numeric(categorical)){
      cat_index <- categorical
    }
  }else{
    cat_index = c()
  }
  
  
  report_list <- list()
  for(i in 1:length(constraints)){  
    if(any_or_not[i] == 0){
      if(constraints_list[i] == "m"){
        report_i <- abs(mumt %*% as.matrix(x[,i]))
      }else if(constraints_list[i] == "s"){
        report_i <- abs((pmt %*% as.matrix(x[,i])))
      }
      if(!i %in% cat_index){
        report_ilist <- quantile(report_i, prob = c(0, 0.25, 0.5, 0.75, 1))
        names(report_ilist) <- c("Minimum", "25th Pctl", "Median", "75th Pctl", "Maximum")
      }else{
        report_ilist <- as.data.frame(table(report_i))
        colnames(report_ilist) <- c("difference", "frequency")
      }
      report_list[[i]] <- report_ilist
    } 
  }
  names(report_list)  <- colnames(x)

  R_result <- dim(qmt)[1]
  

  summary_constraints <- as.data.frame(cbind(S, R, R_result, paste0(round(R_result/R, 4) * 100, "%")))
  colnames(summary_constraints) <- c("overall allocations", "checked allocations", "accepted allocations", "overall % acceptable")
  rw <- sample(R_result, 1)

  inter <- qmt[rw, ]

  }else{

    sim <- 1
    # indicate simulation
    R_over = choose(nsub, ntrt)
    S <- size # randomization sample size

    qmt <- matrix(NA, S, nsub)


    for(s in 1:S){

      trt <- sample(1:nsub, ntrt)
            
      qmt[s, trt] <- 1 # for final
      qmt[s, -trt] <- 0 # for final 
    }


    qmt <- unique(qmt)
    R <- dim(qmt)[1]
  
    pmt <-  omt <- mumt <- omumt <- matrix(NA, R, nsub)

    for(s in 1:R){


      trt <- which(qmt[s, ] == 1)

      omt[s, ] <- 1/2 # for mean arm total
      omumt[s, ] <- 1/nsub  # for overall mean
      
      pmt[s, trt] <- 1 # for arm total
      pmt[s, -trt] <- -1 # for arm total
      
      mumt[s, trt]  <- 1/ntrt # for arm mean
      mumt[s, -trt] <- -1/(nsub - ntrt) # for arm mean
      

    }

  any_or_not <- 1 * (constraints == "any")
  constraints_list <- substring(constraints, 1, 1) # "m" or "s"
  fraction_or_not <- substring(constraints, 2, 2)
  constraints_value <- rep(NA, length(constraints))
  constraints_value[fraction_or_not == "f"] = as.numeric(substring(constraints, 3)[fraction_or_not == "f"])  
  constraints_value[fraction_or_not != "f" & fraction_or_not != "n"] = as.numeric(substring(constraints, 2)[fraction_or_not != "f" & fraction_or_not != "n"])

  use_fraction_or_not = 1 * (fraction_or_not == "f")

  sat_num = R
  for(i in 1:length(constraints)){  
    if(any_or_not[i] == 0){
      if(constraints_list[i] == "m"){
        BL_means <- abs(mumt %*% as.matrix(x[,i]))
        if(use_fraction_or_not[i] == 1){ 
          BL_omean <- omumt %*% as.matrix(x[,i])
          sat_list <- c(1 * (BL_means <= constraints_value[i] * BL_omean))
        }else{
          sat_list <- c(1 * (BL_means <= constraints_value[i]))
        }
      }else if(constraints_list[i] == "s"){
        BL_sums <- abs((pmt %*% as.matrix(x[,i])))
        if(use_fraction_or_not[i] == 1){
          BL_osum <- omt %*% as.matrix(x[,i])
          sat_list <- c(1 * (BL_sums <= constraints_value[i] * BL_osum))
        }else{
          sat_list <- c(1 * (BL_sums <= constraints_value[i]))
        }
      }
    }
    sat_num = sum(sat_list)
    if(sat_num > 0){
      omt <- omt[which(sat_list == 1), ]
      omumt <- omumt[which(sat_list == 1), ]
      pmt <- pmt[which(sat_list == 1), ]
      mumt <- mumt[which(sat_list == 1), ]
      qmt <- qmt[which(sat_list == 1), ]
    }else{
      stop("Error: there is not any scheme that satifies the constraints.")
    }   
  }

  if(!is.null(categorical)){
      if(is.character(categorical)){
        cat_index <- which(colnames(x) %in% categorical)
      }else if(is.numeric(categorical)){
        cat_index <- categorical
      }
    }else{
      cat_index = c()
    }

  report_list <- list()
  for(i in 1:length(constraints)){  
    if(any_or_not[i] == 0){
      if(constraints_list[i] == "m"){
        report_i <- abs(mumt %*% as.matrix(x[,i]))
      }else if(constraints_list[i] == "s"){
        report_i <- abs((pmt %*% as.matrix(x[,i])))
      }
      if(!i %in% cat_index){
        report_ilist <- quantile(report_i, prob = c(0, 0.25, 0.5, 0.75, 1))
        names(report_ilist) <- c("Minimum", "25th Pctl", "Median", "75th Pctl", "Maximum")
      }else{
        report_ilist <- as.data.frame(table(report_i))
        colnames(report_ilist) <- c("difference", "frequency")
      }
        report_list[[i]] <- report_ilist
      } 
  }
  names(report_list)  <- colnames(x)

  R_result <- dim(qmt)[1]
  
  summary_constraints <- as.data.frame(cbind(R_over, R, R_result, paste0(round(R_result/R, 4) * 100, "%")))
  colnames(summary_constraints) <- c("overall allocations", "checked allocations", "accepted allocations", "overall % acceptable")
  rw <- sample(R_result, 1)

  inter <- qmt[rw, ]
   
  }

  if (!is.null(savedata)){

    SchemeChosen <- rep(0, dim(qmt)[1])

    SchemeChosen[rw] <- 1

    pmt <- cbind(SchemeChosen, qmt)

    write.csv(qmt, file = savedata, row.names=FALSE)
  }

  coin_matrix <- coin_descri <- al_clusters <- alno_clusters  <- hi_clusters <- lo_clusters <- NULL
  if(check_validity){
      
      n_pair <- t(combn(nsub, 2))      # all the schemes

      same_arm_count <- same_arm_frac <- diff_arm_count <- diff_arm_frac <- rep(NA, dim(n_pair)[1])
      
      for(j in 1:(nsub-1)){
        for(k in (j+1):nsub){
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
      if(sum(same_arm_frac == 1) > 0){
        alto_index <- which(same_arm_frac == 1)
        alto_clt_pair <- c()
        for(t in 1:length(alto_index)){
          clt_index <- alto_index[t]
          alto_clts <- paste0(first_cluster[clt_index], " and ", second_cluster[clt_index])
          alto_clt_pair <- c(alto_clt_pair, alto_clts)
        }

        alto_all_pt <- rep("100.0%", length(alto_index))
        
        al_clusters <- as.data.frame(cbind(alto_clt_pair, alto_all_pt))
        colnames(al_clusters) <- c("cluter pair", "% allocs in the same arm")
      }

      # Always not together
      if(sum(same_arm_frac == 0) > 0){
        alnoto_index <- which(same_arm_frac == 0)
        alnoto_clt_pair <- c()
        for(t in 1:length(alnoto_index)){
          clt_index <- alnoto_index[t]
          alnoto_clts <- paste0(first_cluster[clt_index], " and ", second_cluster[clt_index])
          alnoto_clt_pair <- c(alnoto_clt_pair, alnoto_clts)
        }

        alnoto_all_pt <- rep("0.0%", length(alnoto_index))
        
        alno_clusters <- as.data.frame(cbind(alnoto_clt_pair, alnoto_all_pt))
        colnames(alno_clusters) <- c("cluter pair", "% allocs in the same arm")
      }
      
      
      # user specify upper bound
      if(sum(same_arm_frac >= samearmhi) > 0){
        hi_index <- which(same_arm_frac >= samearmhi)
        hi_pair <- c()
        for(t in 1:length(hi_index)){
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
      if(sum(same_arm_frac <= samearmlo) > 0){
        lo_index <- which(same_arm_frac <= samearmlo)
        lo_pair <- c()
        for(t in 1:length(lo_index)){
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
  Allocation <- cbind(id, inter)

  colnames(Allocation) <- c("id", "allocation")

  assignment_message <- paste("You have indicated that you want to assign", ntrt, "clusters to treatment", "and", n - ntrt, "to control")


  # indicate a enumeration process or a simulation process with the detailed number of schemes
    if(sim == 1){

    scheme_message <- paste("Simulating", S, "schemes with", R, "unique schemes for", ntrt_cluster, "clusters in the treatment arm out of", ntotal_cluster, "clusters in total")

          } else {

    scheme_message <- paste("Enumerating all the", R, "schemes for", ntrt, "clusters in the treatment arm out of", n, "clusters in total")

    }

       data_merge <- data.frame(inter, id, data)
      # data frame including the chosen scheme from BL, cluster id and the covariates

      colnames(data_merge)[1] <- "arm"

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

      Descriptive_statistics <- CreateTableOne( vars = colnames(data_merge)[c(-1, -2)], strata=c("arm"), data=data_merge, test =  FALSE)
       # create the descriptive table to compare the two arms from constrained randomization using BL

      invisible(capture.output(DS <- as.data.frame(print(Descriptive_statistics))))
     # transform the table into a data frame

      colnames(DS)[1:2] <- c("arm = 0", "arm = 1")

  return(list(allocation = Allocation,
             assignment_message = assignment_message,
             scheme_message = scheme_message,
             data_CR = data_merge,
             baseline_table = DS, 
             cluster_coincidence = coin_matrix, 
             cluster_coin_des = coin_descri, 
             clusters_always_pair = al_clusters, 
             clusters_always_not_pair = alno_clusters, 
             clusters_high_pair = hi_clusters, 
             clusters_low_pair = lo_clusters, 
             overall_allocations = summary_constraints, 
             overall_summary = report_list))

}

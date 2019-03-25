## ----start,echo=FALSE,results="hide"-------------------------------------

library(cvcrand)



## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(Dickinson_design[ , 1:6])

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(Dickinson_design[ , 7:11])

## ----cvrall, fig.keep="all", fig.width = 7, fig.height=4-----------------

 Design_result <- cvrall(clustername = Dickinson_design$county,
                  balancemetric = "l2",
                  x = data.frame(Dickinson_design[ , c("location", "inciis",
                      "uptodateonimmunizations", "hispanic", "incomecat")]),
                  ntotal_cluster = 16,
                  ntrt_cluster = 8,
                  categorical = c("location", "incomecat"),
                  savedata = "dickinson_constrained.csv",
                  bhist = TRUE,
                  cutoff = 0.1,
                  seed = 12345)
 


## ----set-options1, echo=FALSE, fig.keep="all", fig.width = 7, fig.height=4------------------------
options(width = 100)
 

## ---- fig.keep="all", fig.width = 7, fig.height=4-------------------------------------------------
 # the balance metric used
 Design_result$balancemetric

 # the allocation scheme from constrained randomization
 Design_result$allocation
 
 # the histogram of the balance score with respect to the balance metric
 Design_result$bscores
 
 # the statement about how many clusters to be randomized to the intervention and the control arms respectively
 Design_result$assignment_message
 
 # the statement about how to get the whole randomization space to use in constrained randomization
 Design_result$scheme_message
 
 # the statement about the cutoff in the constrained space
 Design_result$cutoff_message
 
 # the statement about the selected scheme from constrained randomization
 Design_result$choice_message
 
 
 # the data frame containing the allocation scheme, the clustername as well as the original data frame of covariates
 Design_result$data_CR
 
 # the descriptive statistics for all the variables by the two arms from the selected scheme
 Design_result$baseline_table

 # the cluster pair descriptive, which is useful for valid randomization check
 Design_result$cluster_coin_des

 # the overall allocation summary
 Design_result$overall_allocations


## ----cvrallst1, fig.keep="all", fig.width = 7, fig.height=4---------------------------------------
# Stratification on location, with constrained randomization on other specified covariates

Design_stratified_result1 <- cvrall(clustername = Dickinson_design$county,
                                     balancemetric = "l2",
                                     x = data.frame(Dickinson_design[ , c("location", "inciis", 
                                                                          "uptodateonimmunizations", 
                                                                          "hispanic", "incomecat")]),
                                     ntotal_cluster = 16,
                                     ntrt_cluster = 8,
                                     categorical = c("location", "incomecat"),
                                     weights = c(1000, 1, 1, 1, 1),
                                     cutoff = 0.1,
                                     seed = 12345) 


## ---- fig.keep="all", fig.width = 7, fig.height=4-------------------------------------------------
Design_stratified_result1$baseline_table


## ----cvrallst2, fig.keep="all", fig.width = 7, fig.height=4---------------------------------------
# An alternative and equivalent way to stratify on location

Design_stratified_result2 <- cvrall(clustername = Dickinson_design$county,
                                     balancemetric = "l2",
                                     x = data.frame(Dickinson_design[ , c("location", "inciis",
                                                                          "uptodateonimmunizations", 
                                                                          "hispanic", "incomecat")]),
                                     ntotal_cluster = 16,
                                     ntrt_cluster = 8,
                                     categorical = c("location", "incomecat"),
                                     stratify = "location",
                                     cutoff = 0.1,
                                     seed = 12345, 
                                     check_validity = TRUE)


## ---- fig.keep="all", fig.width = 7, fig.height=4-------------------------------------------------
Design_stratified_result2$baseline_table


## ----cvrcov, fig.keep="all", fig.width = 7, fig.height=4------------------------------------------

# change the categorical variable of interest to have numeric representation
Dickinson_design_numeric <- Dickinson_design
Dickinson_design_numeric$location = (Dickinson_design$location == "Rural") * 1

Design_cov_result <- cvrcov(clustername = Dickinson_design_numeric$county,
                            x = data.frame(Dickinson_design_numeric[ , c("location", "inciis", 
                                                                          "uptodateonimmunizations", 
                                                                          "hispanic", "income")]),
                            ntotal_cluster = 16,
                            ntrt_cluster = 8,
                            constraints = c("s5", "mf.5", "any", "any", "mf0.4"), 
                            categorical = c("location"),
                            savedata = "dickinson_cov_constrained.csv",
                            seed = 12345, 
                            check_validity = TRUE)
 


## ----set-options2, echo=FALSE, fig.keep="all", fig.width = 7, fig.height=4------------------------
options(width = 100)
 

## ---- fig.keep="all", fig.width = 7, fig.height=4-------------------------------------------------


 # the allocation scheme from constrained randomization
 Design_cov_result$allocation
 

 # the statement about how many clusters to be randomized to the intervention and the control arms respectively
 Design_cov_result$assignment_message
 
 # the statement about how to get the whole randomization space to use in constrained randomization
 Design_cov_result$scheme_message
 

 # the data frame containing the allocation scheme, the clustername as well as the original data frame of covariates
 Design_cov_result$data_CR
 
 # the descriptive statistics for all the variables by the two arms from the selected scheme
 Design_cov_result$baseline_table

# the cluster pair descriptive, which is useful for valid randomization check
Design_cov_result$cluster_coin_des

# the overall allocation summary
Design_cov_result$overall_allocations


## ---- echo=FALSE, results='asis'------------------------------------------------------------------
knitr::kable(head(Dickinson_outcome, 10))

## ----cptest, fig.keep="all", fig.width = 7, fig.height=4------------------------------------------
 Analysis_result <- cptest(outcome = Dickinson_outcome$outcome,
                           clustername = Dickinson_outcome$county,
                           z = data.frame(Dickinson_outcome[ , c("location", "inciis",
                               "uptodateonimmunizations", "hispanic", "incomecat")]), 
                            cspacedatname = system.file("dickinson_constrained.csv", package = "cvcrand"),                                 
                           outcometype = "binary",                                                      
                           categorical = c("location","incomecat"))



## ----cptestre, fig.keep="all", fig.width = 7, fig.height=4----------------------------------------
 Analysis_result 

## ----info, results='markup', echo=FALSE-----------------------------------------------------------
sessionInfo()


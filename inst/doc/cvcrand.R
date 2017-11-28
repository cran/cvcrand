## ----setup,echo=FALSE,results="hide"-------------------------------------

library(cvcrand)



## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(Dickinson_design)

## ----cvcrand, fig.keep="all", fig.width = 7, fig.height=4----------------

 Design_result<-cvcrand(id = Dickinson_design$county,
                       metric = "L2",
                       x = data.frame(Dickinson_design[, -1]),
                       n = 16,
                       ntrt = 8,
                       categorical = c("location","incomecat"),
                       savedata = "dickinson_constrained.csv",
                       savebscores = "dickinson_bscores.csv",
                       cutoff = 0.1,
                       seed = 12345)
 
 



## ----set-options, echo=FALSE, fig.keep="all", fig.width = 7, fig.height=4-------------------------
options(width = 100)
 

## ---- fig.keep="all", fig.width = 7, fig.height=4-------------------------------------------------
 # the metric the user specified
 Design_result$metric

 # the selected scheme 
 Design_result$Allocation
 
 # the cutoff balance score, the balance score of the selected   scheme, as well as the histogram of the balance scores of the whole randomization space
 Design_result$Bscores
 
 # the statement about how many clusters to randomize to the intervention and the control arms respectively
 Design_result$assignment_message
 
 # the statement about how to get the whole space of schemes
 Design_result$scheme_message
 
 # the statement about the benchmark of the constrained space
 Design_result$BLcut_message
 
 # the statement about the selected scheme from constrained randomization
 Design_result$BLchoice_message
 
 
 # the matrix containing allocation scheme, the id as well as the original covariates' matrix
 Design_result$data_BL
 
 # the descriptive statistics for all the variables in the original covariates' matrix in the two arms from constrained randomization
 Design_result$BL_result



## ---- echo=FALSE, results='asis'------------------------------------------------------------------
knitr::kable(head(Dickinson_outcome, 10))

## ----cptest, fig.keep="all", fig.width = 7, fig.height=4------------------------------------------
 Analysis_result <- cptest(outcome = Dickinson_outcome$outcome,
                           id = Dickinson_outcome$county,
                           x = data.frame(Dickinson_outcome[, c(-1, -7)]),
                           cspacedatname = system.file("dickinson_constrained.csv", package="cvcrand"),
                           outcometype = "binary",
                           categorical = c("location","incomecat"))


## ----cptestre, fig.keep="all", fig.width = 7, fig.height=4----------------------------------------
 Analysis_result 

## ----info, results='markup', echo=FALSE-----------------------------------------------------------
sessionInfo()


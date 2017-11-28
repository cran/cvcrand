[![Build Status](https://travis-ci.com/hengshiyu/cvcrand.svg?token=ZtDwjW3Z5FqDgnzwZCn6&branch=master)](https://travis-ci.com/hengshiyu/cvcrand)

# cvcrand R-package
Hengshi Yu, John A. Gallis, Fan Li and Elizabeth L. Turner

**Maintainer**: Hengshi Yu <hengshi@umich.edu>

## Introduction
cvcrand is an R package for the design and analysis of cluster randomized trials (CRTs). 

A cluster is the unit of randomization for a cluster randomized trial. Thus, when the number of clusters is small, there might be some baseline imbalance from the randomization between the arms. Constrained randomization constrained the randomization space. Given the baseline values of some cluster-level covariates, users can perform a constrained randomization on the clusters into two arms, with an optional input of user-defined weights on the covariates. 

At the end of the study, the individual outcome is collected. The `cvcrand` package then performs clustered permutation test on either continuous outcome or binary outcome adjusted for some individual-level covariates, producing p-value of the intervention effect.

## Functions and references
The cvcrand package constains two main functions. In the design of CRTs with two arms, users can use the `cvcrand()` function to perform constrained randomization. And for the analysis part, user will use the `cptest()` function for clustered permutation test. 

1. cvcrand function: constrained randomization for two-arm cluster randomized trials
    * Raab, G.M. and Butcher, I., 2001. Balance in cluster randomized trials. Statistics in medicine, 20(3), pp.351-365.

2. cptest function: clustered permutation test for two-arm cluster randomized trial
    * Gail, M.H., Mark, S.D., Carroll, R.J., Green, S.B. and Pee, D., 1996. On design considerations and randomization‐based inference for community intervention trials. Statistics in medicine, 15(11), pp.1069-1092.
    * Li, F., Lokhnygina, Y., Murray, D.M., Heagerty, P.J. and DeLong, E.R., 2016. An evaluation of constrained randomization for the design and analysis of group‐randomized trials. Statistics in medicine, 35(10), pp.1565-1579.

### `cvcrand()` example: constrained randomization

The balance score for constrained randomization in the program is developed from Raab and Butcher (2001).  

Study 1 in Dickinson et al (2015) is about two approaches (interventions) for increasing the "up-to-date" immunization rate in 19- to 35-month-old children. They planned to randomize 16 counties in Colorado 1:1 to either a population-based approach or a practice-based approach. There are several county-level variables. The program will randomize on a subset of these variables. The continuous variable of average income is categorized to illustrate the use of the `cvcrand()` on multi-category variables. And the percentage in CIIS variable is truncated at 100%.

For the constrained randomization, we used the `cvcrand()` function to randomize 8 out of the 16 counties into the practice-based. For the definition of the whole randomization space, if the total number of all possible schemes is smaller than `50,000`, we enumerate all the schemes as the whole randomization space. Otherwise, we simulate `50,000` schemes and choose the unique shemes among them as the whole randomization space. We calculate the balance scores of `"L2"` metric on three continuous covariates as well as two categorical covariates of location and income category. Location has `"Rural"` and `"Urban"`. The level of `"Rural"` was then dropped in `cvcrand()`. As income category has three levels of `"low"`, `"med"`, and `"high"`,  the level of `"high"` was dropped to create dummy variables according to the alphanumerical order as well. Then we constrained the randomization space to the schemes with `"L2"` balance scores less than the `0.1` quantile of that in the whole randomization space. Finally, a randomization scheme is sampled from the constrained space. 

We saved the constrained randomization space in a CSV file in `"dickinson_constrained.csv"`, the first column of which is an indicator variable of the finally selected scheme (`1`) or not (`0`). We also saved the balance scores of the whole randomization space in a CSV file in `"dickinson_bscores.csv"`, and output a histogram displaying the distribution of all balance scores with a red line indicating our selected cutoff (the `0.1` quantile).  



```r

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


```

### `cptest()` example: Clustered Permutation Test

At the end of cluster randomized trials, individual outcomes are collected. Permutation test based on Gail et al (1996) and Li et al (2016) is then applied to the continuous or binary outcome with some individual-level covariates. 

Suppose that the researchers were able to assess 300 children in each cluster of the study 1 in [@dickinson2015pragmatic], and the cluster randomized trial is processed with the selected randomization scheme from the example above of the `cvcrand()` function. We expanded the values of the cluster-level covariates on the covariates' values of the individuals, according to which cluster they belong to. The correlated individual outcome of up-to-date on immunizations (`1`) or not (`0`) is then simulated using a generalized linear mixed model (GLMM) to induce correlation by include a random effect. The intracluster correlation (ICC) was set to be 0.01, using the latent response definition provided in [@eldridge2009intra]. This is a reasonable value of the ICC the population health studies [@hannan1994parameters]. We simulated one data set, with the outcome data dependent on the county-level covariates used in the constrained randomization design and a positive treatment effect so that the practice-based intervention increases up-to-date immunization rates more than the community-based intervention. For each individual child, the outcome is equal to `1` if he or she is up-to-date on immunizations and `0` otherwise.with differential sampling by intervention arm. 

We used the `cptest()` function to process the clustered permutation test on the binary outcome of the status of up-to-date on immunizations. We input the file about the constrained space with the first column indicating the final scheme. The permutation test is on the continuous covariates of `"inciis"`, `"uptodateonimmunizations"`, `"hispanic"`, as well as categorical variables of `"location"` and `"incomecat"`. Location has `"Rural"` and `"Urban"`. The level of `"Rural"` was then dropped in `cptest()`. As income category has three levels of `"low"`, `"med"`, and `"high"`,  the level of `"high"` was dropped to create dummy variables according to the alphanumerical order as well.

```r
 Analysis_result <- cptest(outcome = Dickinson_outcome$outcome,
                           id = Dickinson_outcome$county,
                           x = data.frame(Dickinson_outcome[, c(-1, -7)]),
                           cspacedatname = system.file("dickinson_constrained.csv", package="cvcrand"),
                           outcometype = "binary",
                           categorical = c("location","incomecat"))
```

## Installation

This is the first time the `cvcrand` R package is submitted to CRAN.

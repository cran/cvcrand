# cvcrand 0.0.4

## Major changes

1. changed the default `size` of `cvrcov()` function to be `50,000` from `100,000` and to simulate `size` schemes and remove duplicates like in `cvrall()` function rather than previously directly `size` unique schemes. 
2, revised both `cvrall()` and `cvrcov()` for keeping `cluster_coin_des` to the third decimal. 

3, revised the clustered permutation test in `cptest()` function to use the residuals in the `response` level rather than in the `logit` level. 

4, revised `cvrcov()` function for its output. We put `"arm"` as the first column of `data_CR` meaning the allocated arm and also changed column names of `baseline_table` to be `arm = 0` and `arm = 1`
